#!/bin/bash

set -e
# This wrapper script dynamically creates volume mounts for any provided host path.

# --- Correctly find and create the output directory ---
HOST_OUTPUT_PATH=""
ARGS=("$@") 
for arg in "${ARGS[@]}"; do
    if [[ "$arg" == "--help" || "$arg" == "-h" ]]; then
        echo "Forwarding help request to the SmedAnno pipeline..."
        # Create a dummy .env file as docker-compose requires it,
        # but the path doesn't matter for a help request.
        echo "OUTPUT_PATH=/tmp/smedanno_dummy_output" > .env
        docker-compose run --rm smedanno --help
        rm .env # Clean up
        exit 0
    fi
done
for i in "${!ARGS[@]}"; do
    if [[ "${ARGS[$i]}" == "--outputDir" ]]; then
        j=$((i + 1))
        if [[ $j -lt ${#ARGS[@]} ]]; then
            # Just store the path as a string for now
            HOST_OUTPUT_PATH="${ARGS[$j]}"
            break 
        fi
    elif [[ "${ARGS[$i]}" == --outputDir=* ]]; then
        HOST_OUTPUT_PATH="${ARGS[$i]#*=}"
        break
    fi
done

if [ -z "$HOST_OUTPUT_PATH" ]; then
    echo -e "\033[0;31mError: You must provide the --outputDir argument.\033[0m"
    exit 1
fi

# Create the directory
#mkdir -p "$HOST_OUTPUT_PATH"
#HOST_OUTPUT_PATH=$(realpath "$HOST_OUTPUT_PATH")

# Handle the output directory creation and validation
if [ -e "$HOST_OUTPUT_PATH" ] && [ ! -d "$HOST_OUTPUT_PATH" ]; then
    # Path exists but it's a file, not a directory. This is a fatal error.
    echo -e "\033[0;31mError: Output path '$HOST_OUTPUT_PATH' exists but is a file. Please remove it or choose a different path.\033[0m"
    exit 1
fi

# Attempt to create the directory. The -p flag ensures no error if it already exists.
if ! mkdir -p "$HOST_OUTPUT_PATH"; then
    echo -e "\033[0;31mError: Failed to create output directory '$HOST_OUTPUT_PATH'.\033[0m"
    echo "Please check write permissions for the parent directory."
    exit 1
fi

# Now that we are certain the path exists and is a directory, resolve it.
HOST_OUTPUT_PATH=$(realpath "$HOST_OUTPUT_PATH")
if [ $? -ne 0 ]; then
     echo -e "\033[0;31mError: Could not resolve the absolute path for the output directory. Please check its validity.\033[0m"
     exit 1
fi

# Check if the directory is writable
if [ ! -w "$HOST_OUTPUT_PATH" ]; then
    echo -e "\033[1;33mWarning: Output directory is not writable by the current user ($(whoami)).\033[0m"
    echo "Attempting to adjust permissions on $HOST_OUTPUT_PATH to allow container access (chmod o+w)..."
    if ! chmod o+w "$HOST_OUTPUT_PATH"; then
        echo -e "\033[0;31mError: Failed to set write permissions for the output directory.\033[0m"
        echo "Please run 'sudo chmod o+w $HOST_OUTPUT_PATH' and try again."
        exit 1
    fi
fi

# --- Initialize variables for dynamic mounting ---
VOLUME_ARGS=""
MOUNT_COUNTER=0
declare -A DIR_MAP
final_container_args=()

# --- Main argument processing loop ---
i=0
while [ $i -lt ${#ARGS[@]} ]; do
    arg="${ARGS[$i]}"
    
    case "$arg" in
        --genomeRef|--genomeGTF|--rrnaRef|--finalGTF|--alignDirShort|--alignDirMix|--dataDirShort|--dataDirMix|--SR_RB_gtf_dir|--SR_DN_gtf_dir|--MR_RB_gtf_dir|--MR_DN_gtf_dir|--genomeDir)
            final_container_args+=("$arg")
            i=$((i + 1))
            
            POTENTIAL_HOST_PATH="${ARGS[$i]}"
            if [ ! -e "$POTENTIAL_HOST_PATH" ]; then
                echo -e "\033[0;31mError: Input path not found for argument $arg: $POTENTIAL_HOST_PATH\033[0m"
                exit 1
            fi

            HOST_PATH=$(realpath "${POTENTIAL_HOST_PATH}")
            
            if [ -d "$HOST_PATH" ]; then
                HOST_DIR="$HOST_PATH"; FILENAME=""
            else
                HOST_DIR=$(dirname "$HOST_PATH"); FILENAME=$(basename "$HOST_PATH")
            fi

            if [[ -z "${DIR_MAP[$HOST_DIR]}" ]]; then
                MOUNT_COUNTER=$((MOUNT_COUNTER + 1))
                CONTAINER_MOUNT_POINT="/mount${MOUNT_COUNTER}"
                DIR_MAP[$HOST_DIR]=$CONTAINER_MOUNT_POINT
                VOLUME_ARGS+=" -v ${HOST_DIR}:${CONTAINER_MOUNT_POINT}:ro"
            fi
            
            CONTAINER_DIR="${DIR_MAP[$HOST_DIR]}"
            CONTAINER_PATH=$([ -n "$FILENAME" ] && echo "${CONTAINER_DIR}/${FILENAME}" || echo "${CONTAINER_DIR}")
            final_container_args+=("$CONTAINER_PATH")
            ;;
        
        --outputDir)
            final_container_args+=("$arg"); final_container_args+=("/output"); i=$((i + 1))
            ;;

        *)
            final_container_args+=("$arg")
            if [[ $i+1 -lt ${#ARGS[@]} ]] && ! [[ ${ARGS[$i+1]} =~ ^-- ]]; then
                i=$((i + 1)); final_container_args+=("${ARGS[$i]}")
            fi
            ;;
    esac
    i=$((i + 1))
done

# Create .env file
echo "Creating dynamic .env file..."
echo "OUTPUT_PATH=${HOST_OUTPUT_PATH}" > .env

# Execute Docker Compose
echo "Executing pipeline via docker-compose with dynamic mounts..."
docker-compose run --rm ${VOLUME_ARGS} smedanno "${final_container_args[@]}"

echo "Wrapper script finished."

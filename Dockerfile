# Start from pre-built base image
FROM smedanno-base

# Set the working directory
WORKDIR /smedanno

# Copy scripts and set permissions
COPY --chmod=755 smedanno.sh functional_annotation.R tagXSstrandedData.awk ./
#RUN chmod +x smedanno.sh

# Define the entrypoint to run inside the smedanno environment
ENTRYPOINT ["conda", "run", "--live-stream", "-n", "smedanno", "/smedanno/smedanno.sh"]

# !/bin/bash
echo "Stopping and removing container"
docker stop thesis
docker rm thesis
echo "Done!"

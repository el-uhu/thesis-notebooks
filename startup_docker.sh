# !/bin/bash
echo "Starting container and entering shell"
docker run -d --name thesis -p 8888:8888 eluhu/thesis-nb
docker ps -l
echo "Jupyter notebook running under http://localhost:8888"
docker exec -it thesis /bin/bash

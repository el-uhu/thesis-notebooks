# About this repository

This repository contains interactive notebooks that document the models I used in the work covered in my DPhil thesis, and make them fully interactive.
Thus, parameters can be changed very easily, and simulations can be re-run.

The code I developed has numerous dependencies, such as the programming language `julia`,  the dynamical systems modelling programm `xpp` and various packages to enable plotting and fitting.
It can be a hassle to install these on different systems, as dependencies can clash, or crucial dependencies can be omitted from the documentation.

To circumvent this problem, I set up a `docker` container. Docker containers can be thought of as portable self-contained computing environments. They work like slimmed-down virtual operating systems, that can be run like a program on your computer. That way, none of the packages required to replicate my simulations will touch your operating system, but things will work smoothly, as everything my software needs is installed within the container.
If you are interested in having a look at the inside of this custom docker container, go to: http://github.com/el-uhu/xyz. There you will find a so-called `Dockerfile`. Much like a recipe, this contains the instructions on how your computer builds the container in a step by step manner.
These steps include cloning of this repository. Therefore, these instructions, as well as the interactive notebooks will be immediately available to you, when you run the container.

If you want to check out the notebooks and don't need them to be fully interactive, just navigate to the subfolders and click on the files, and you will see a preview.

# Installing docker

## Mac OS
Follow the instructions on https://docs.docker.com/docker-for-mac/

Requirements:
- Mac must be a 2010 or newer model, with Intel’s hardware support for memory management unit (MMU) virtualization; i.e., Extended Page Tables (EPT)
- OS X 10.10.3 Yosemite or newer
- At least 4GB of RAM
- VirtualBox prior to version 4.3.30 must NOT be installed (it is incompatible with Docker for Mac). Docker for Mac will error out on install in this case. Uninstall the older version of VirtualBox and re-try the install.

## Windows
Follow the instructions on https://docs.docker.com/docker-for-windows/

## Ubuntu Linux (12.04 or higher)
Follow the instructions on https://docs.docker.com/engine/installation/linux/ubuntulinux/

## Other Linux Distributions
Look for suitable documentation here: https://docs.docker.com/engine/installation/

# Installing, Running, and Stopping  the Custom Docker Container


## Install
Get the docker container:
```
docker pull eluhu/jupyter-julia
```
It contains all the software packages, code and data needed to run simulations interactively.

## Run
Once it is downloaded, run it using:

```
docker run -d --name thesis -p 8888:8888 eluhu/jupyter-julia jupyter notebook
```

This command will start up an interactive notebook in the background. It will have the name thesis.

To check, whether it is running correctly, type:
```
docker ps -a
```
This will list all running container instances.
The output should contain the following line:
```
CONTAINER ID        IMAGE                 COMMAND              CREATED             STATUS                         PORTS                    NAMES
e1edaa91217f        eluhu/jupyter-julia   "jupyter notebook"   3 seconds ago       Up 2 seconds                   0.0.0.0:8888->8888/tcp   thesis
```
Importantly, `STATUS` should show `Up`, and `PORTS` should display  `0.0.0.0:8888->8888/tcp`.

To use the notebooks, open http://localhost:8888 in a browser.
You will see folders, corresponding to the chapters of my thesis on results of my research.

## Stop
To stop the container named `thesis` type:
```
docker stop thesis
```

You can use `docker ps -a` to make sure that the `STATUS` of `thesis` has changed to `Exited (0) 1 seconds ago`

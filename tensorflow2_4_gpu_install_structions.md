# Install Tensorflow 2.4

[Github pages version](https://markjay4k.github.io/Install-Tensorflow/)

We're going for [Tensorflow 2.4.1 with GPU support](https://www.tensorflow.org/install/gpu) on an Ubuntu 20.04 system.
- Sorry, but we won't use the [Tensorflow Docker Image with GPU Support](https://www.tensorflow.org/install/docker).

My hardware specs are:
- GPU: GTX 1070
- OS: Ubuntu 20.04

## Step 1 - Update Nvidia GPU Driver

We will install using the PPA Repository. This will install the latest driver. Run the following command to install the `graphics-driver/ppa` repository to your system

```shell
sudo add-apt-repository ppa:graphics-drivers/ppa
```
Then install the latest driver. Though we call out `440`, it should install the latest (`450` in my case)

```shell
sudo apt install nvidia-driver-440
```
You need to reboot to complete the driver installation.

```shell
sudo reboot
```
Check that the driver is Install

```shell
nvidia-smi
```
## Step 2 - Install CUDA Toolkit

Download [CUDA Toolkit 11.0](https://developer.nvidia.com/cuda-11.0-download-archive)
- Select Linux, x86_64, Ubuntu, 20.04, deb (local)

Run the commands provided by Nvidia (also shown below).
- NOTE: The toolkit about 2.4GB

```shell
wget https://developer.download.nvidia.com/compute/cuda/repos/ubuntu2004/x86_64/cuda-ubuntu2004.pin
sudo mv cuda-ubuntu2004.pin /etc/apt/preferences.d/cuda-repository-pin-600
wget http://developer.download.nvidia.com/compute/cuda/11.0.2/local_installers/cuda-repo-ubuntu2004-11-0-local_11.0.2-450.51.05-1_amd64.deb
sudo dpkg -i cuda-repo-ubuntu2004-11-0-local_11.0.2-450.51.05-1_amd64.deb
sudo apt-key add /var/cuda-repo-ubuntu2004-11-0-local/7fa2af80.pub
sudo apt update
sudo apt install cuda
```

If you get an error running the last command to install, you may have conflicting CUDA packages in your `/etc/apt/preferences.d` directory. Go in and remove the unwanted ones if this occurs.

## Step 4 - Install CUDNN 8.0.4

You need an Nvidia account to download cuDNN, so sign up for one (don't worry, it's free). Download cuDNN v8.0.4 Library for Linux (x86_64) from the [cuDNN Archives](https://developer.nvidia.com/rdp/cudnn-archive#a-collapse804-110).
Once, downloaded, we need to unzip the tar file (`tar`), copy some files to our CUDA Toolkit folders (`cp`), and change their mode (`chmod`). To do this, run the following commands:

```shell
tar -xzvf cudnn-11.0-linux-x64-v8.0.4.30.tgz
sudo cp cuda/include/cudnn.h /usr/local/cuda/include
sudo cp cuda/lib64/libcudnn* /usr/local/cuda/lib64
sudo chmod a+r /usr/local/cuda/include/cudnn.h /usr/local/cuda/lib64/libcudnn*
```

Next, we need to update our `.bashrc` (or `.zshrc` if you're using `zsh`) by adding the path to a few cuda folders
to `LD_LIBRARY_PATH`. Add the following lines to the end of our `.bashrc` file (or `.zshrc` file)

```shell
export LD_LIBRARY_PATH=/usr/local/cuda-11.0/lib64:$LD_LIBRARY_PATH
export LD_LIBRARY_PATH=/usr/local/cuda/lib64:$LD_LIBRARY_PATH
```

## Step 5 - PIP Install

We will be installing Tensorflow using pip3, which is not installed by default with Ubuntu. You can installed using apt

```shell
sudo apt update
sudo apt install python3-pip
```

We will be using virutal-env, which is also not installed by default with python3 on
Ubuntu. If needed, install `venv` with apt.

```shell
sudo apt install python3-venv
```

We will create a new Virtual Environment with `venv` and install Tensorflow
there. First, let's make a new virtual environment.
- NOTE: the environment will be the same version of python that your system
  uses. In my case, this is python-3.8

I like to name my virtual environments `venv` and place them in a directory named after the project name. For example, I create a folder called `tensorflow` and then create an environment called `venv` within the `tensorflow` directory.

```shell
mkdir tensorflow
cd tensorflow
python3 -m venv venv
```

Nex, activate the environment

```shell
source venv/bin/activate
```

Now install tensorflow with pip.

```shell
pip install tensorflow
```

## Step 6 - Test It!

We should be good to go at this point. With the environment still activated, start python, import tensorflow, and run the command below to check for GPUs on your system.

```python
import tensorflow as tf
tf.config.list_physical_devices('GPU')
```

you should see a bunch of output lines and finally a list of your devices. In my case, I get

```python
[PhysicalDevice(name='/physical_device:GPU:0', device_type='GPU')]
```

If you get a device listed, you're good! If you get an empty list, something went wrong. Most likely there is an issue with the CUDNN files not being copied over correctly, or your `LD_LIBRARY_PATH` is incorrect.
that's it! Enjoy!!

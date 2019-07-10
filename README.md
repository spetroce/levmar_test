# Introduction
This is the repository for the Wash Tank firmware. You may import this as a project into Atollic TrueSTUDIO or build it using CMake. All MCU configurations should be done through CubeMX.

# Pull the repo with submodules
```bash
git clone --recurse-submodules git@github.com:makerbot/wash_tank.git
```
*or*
```bash
git clone git@github.com:makerbot/wash_tank.git
cd wash_tank #ie. change to root repo directory
git submodule update --init --recursive
```
If you do not have a github account setup with an SSH key, which is likely the case for EE and ME team members, replace the URL with https://github.com/makerbot/wash_tank.git. So, the shorthand clone instruction is:
```bash
git clone --recurse-submodules https://github.com/makerbot/wash_tank.git
```

# Build Prerequisites
Prerequisite packages to install to perform any building:

Ubuntu:
```bash
sudo apt-get install cmake build-essential git-lfs
```
Fedora:
```bash
sudo dnf install make automake cmake gcc gcc-c++ git git-lfs libtool
```
If you want to build the firmware, you need to have the arm-none-eabi-gcc, arm-none-eabi-g++, arm-none-eabi-copy, arm-none-eabi-size command line tools.

Ubuntu:
```bash
sudo apt-get install gcc-arm-none-eabi
```
Fedora:
```bash
sudo dnf install arm-none-eabi-newlib arm-none-eabi-gcc
```
# Build Firmware and Host Driver

Then, create a build directory:
```bash
cd wash_tank
mkdir build
cd build
```
If you do not wish to build the firmware (ie. only build the host driver controller), issue the following command:
```bash
cmake -DBUILD_FIRMWARE=OFF ..
```
To build both the firmware and the host driver:
```bash
cmake ..
```
Then, call make to build:
```bash
make
```

# Flash Firmware to MCU
You can flash directly to the MCU through an ST-Link debugger.

## Prerequisites
Install the stlink utilities from your package manager or by building/installing the [open source stlink](https://github.com/texane/stlink). 

Ubuntu:

Ubuntu does not have stlink available through the package manager as of writing this. For Ubuntu, you must build/install from source.

Fedora:
```bash
sudo dnf install stlink
```
## Perform Firmware Flash
To flash the firmware, simply connect an ST-Link debugger (e.g. an ST Discovery board) to the MCU via SWD and to your PC via USB. Then run the following command:
```bash
cd wash_tank/build
make install
```

# Run Host Driver (controller)
You can communicate with the microcontroller over UART with the compiled host driver controller called wash_tank_driver. As of now you can set the PWM duty cycle for the motor and read back the current RPM. More features will be added in the future.
```bash
cd wash_tank/build
./modules/driver/wash_tank_driver /dev/ttyUSB0
```
The USB device will typically be /dev/ttyUSB0, but could also exists under any /dev/ttyUSBn. Once the driver has been executed, usage instructions will print to the console.

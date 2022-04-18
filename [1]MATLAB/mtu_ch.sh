#!/bin/sh
sudo ifconfig enp1s0f0 mtu 9000
sudo ifconfig enp1s0f1 mtu 9000

sudo sysctl -w net.core.rmem_max=24862979
sudo sysctl -w net.core.wmem_max=24862979
#sudo ifconfig enp1s0f0 mtu 1500
#sudo ifconfig enp1s0f1 mtu 1500

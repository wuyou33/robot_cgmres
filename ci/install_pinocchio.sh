sudo sh -c "echo 'deb [arch=amd64] http://robotpkg.openrobots.org/packages/debian/pub bionic robotpkg' >> /etc/apt/sources.list.d/robotpkg.list"
curl http://robotpkg.openrobots.org/packages/debian/robotpkg.key | sudo apt-key add -
sudo apt update
sudo apt install robotpkg-py27-pinocchio
echo export PATH=/opt/openrobots/bin:$PATH >> ~/.bashrc
echo export PKG_CONFIG_PATH=/opt/openrobots/lib/pkgconfig:$PKG_CONFIG_PATH >> ~/.bashrc
echo export LD_LIBRARY_PATH=/opt/openrobots/lib:$LD_LIBRARY_PATH >> ~/.bashrc
echo export PYTHONPATH=/opt/openrobots/lib/python2.7/site-packages:$PYTHONPATH >> ~/.bashrc
source ~/.bashrc
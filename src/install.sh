python setup.py sdist --dist-dir='../dist'
cd ../dist
tar zxvf *.tar.gz
cd PyCosmology*
sudo python setup.py install
cd ..
sudo rm -rf PyCosmol*/*
rmdir PyCosmology*
cd ../src

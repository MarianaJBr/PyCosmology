python setup.py sdist --dist-dir='../../PyCosmology_dist'
cd ../../PyCosmology_dist
tar zxvf *.tar.gz
rm *.tar.gz
cd PyCosmology*
sudo python setup.py install
cd ../../PyCosmology/src
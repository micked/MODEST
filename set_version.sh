#Add new version
python2 setup.py --version

#Amend last commit
git add mage_tool/version.py
git commit --amend -C HEAD

#!/bin/bash

# Authors:
#   Wenhan TANG - 08/2021
#   ...

Me=$( readlink -m $( type -p $0 ))
MyDir=`dirname $Me`
cd $MyDir

echo "#!/bin/vim" > README.tree
echo "Bayesian Inversion System v6.0" >> README.tree
echo "Program structure:" >> README.tree

tree >> README.tree

sed -i "/__init__.py/d" README.tree 
sed -i "/__pycache__/d" README.tree
sed -i "/\.pyc/d" README.tree

echo "# Created by tree.sh" >> README.tree
echo "" >> README.tree

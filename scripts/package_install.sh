#!/bin/bash

# You can also install the development version use git clone to
# automatically handle some issues related to dependencies:
git clone https://github.com/mengxu98/inferCSN.git
cd inferCSN
sh scripts/requirements.sh
R CMD INSTALL . --library=/your/lib/path

#!/bin/bash
set -e

echo "🛠️ Update..."
sudo apt update
sudo apt install -y software-properties-common curl wget build-essential
sudo apt install -y python3-pip
pip3 install virtualenv  --break-system-packages



echo "✅ Install gmp、mpfr、qd ..."
sudo apt-get install -y autoconf automake pkg-config
sudo apt-get install -y libtool-bin libgmp-dev
sudo apt-get install -y libmpfr-dev
sudo apt-get install -y libqd-dev

echo "📂 Clone G6K..."
rm -rf g6k
git clone https://github.com/fplll/g6k.git
cd g6k
pip install -r requirements.txt --break-system-packages
PYTHON=python3 ./bootstrap.sh -j 8

alias pydoc 2>/dev/null >/dev/null && unalias pydoc || true
alias python=python3
alias pip="python3 -m pip"


echo "📂  Clone our method..."
rm -rf Modular-hints
git clone https://github.com/Halowooder/Modular-hints.git
cp -r g6k Modular-hints

echo "🎉 Congratulations..."

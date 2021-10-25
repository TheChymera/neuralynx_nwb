# Whole-Project Neuralynx to NWB converter

This is a work in progress script to automate conversion for the “vStr_phase_stim” project by Manish Mohapatra of Dartmouth College.
It is best used from the live `HEAD` as it is subject to active development.

## Installation

Fetch the newest version git and install it using `setuptools` (this will also allow you to use changes you manually make to the repository files):

```
git clone git@github.com:TheChymera/neuralynx_nwb.git
cd neuralynx_nwb
echo "export PATH=\$HOME/.local/bin/:\$PATH" >> ~/.bashrc
source ~/.bashrc
python setup.py develop --user
```

If you are getting a `Permission denied (publickey)` error upon trying to clone, you can either:

* [Add an SSH key](https://help.github.com/articles/adding-a-new-ssh-key-to-your-github-account/) to your GitHub account.
* Pull via the HTTPS link `git clone https://github.com/TheChymera/neuralynx_nwb.git`.

## Usage

The fundamental usage of the data conversion (assuming data is located under the default path) can be initiated via:

```
python -c 'import convert; convert.reposit_data(debug=True)'
```

## Data

The conversion tool currently uses the `vStr_phase_stim` data, which should be located under an eponymous directory inside `~/.local/share/datalad/`, i.e. `~/.local/share/datalad/vStr_phase_stim`.
If the data collecition is too large for the user home directory, or is shared with other users on the machine, it can be located under a different partition (e.g. `/mnt/data/datalad`) and a symlink can be created for `~/.local/share/datalad/` to point to that location.
The example command for these paths would be:

```
ln -s /mnt/data/datalad ~/.local/share/datalad
```

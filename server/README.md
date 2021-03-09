# Backend

## Getting started (Windows)

### Installing Python

Download and install Python from https://www.python.org/downloads/

Python packages will be installed to virtual environments later on. See here for more information: https://docs.python.org/3/tutorial/venv.html

Virtual environments are, in essence, separate Python instances with their own installed packages. It is important to note that once a virtual environment is created, its file path cannot be (easily) changed without breaking the environment. Virtual environments must be activated each time a new console window is opened (assuming one wants to use the virtual environment in that window).

### Cloning the project

Clone the ePotku repository with git (get the correct address from the GitHub repository's Clone button):
```
git clone git@github.com:JYU-IBA/ePotku.git
cd ePotku
```

Checkout the branch you want to use (Master is selected by default):
```
git checkout master
```

### Cloning submodules

Pull, update and initialize submodules:
```
git submodule update --init --recursive
```

### Installing the dependencies

Create (once per project) and activate a virtual environment (every time). In this example, the environment's name is `env`.
```
cd server
py -3 -m venv env
env\Scripts\activate
```

Install the depencies:
```
pip install -r requirements.txt
```

### Installing / compiling the submodules

Refer to Potku's GitHub README for instructions. Explicitly installing Potku itself is not needed, but its accompanying C programs must be compiled. Make sure to follow the right version's instructions. At the time of writing, the correct commit SHA for Potku submodule is eb850fad6c057cb90a0ef65e2b0cac8bb0bf9b0d. The version used for the submodule is visible in GitHub's web interface, in the server folder (here if you are reading this on GitHub).

### Running the application

Run the application (remember to activate your virtual environment if you haven't already):
```
env\Scripts\activate
python app.py
```

Alternative way to run the application:
```
env\Scripts\activate
python -m flask run
```

----

## Getting started (Linux, untested)

The fundamental commands can be the same as on Windows. The biggest difference is in installing and running Python. Virtual environment commands are slightly different too.

The best way to install Python depends on your Linux distribution. There are many possible choices, but using Pyenv to install a specific version is probably a safe option. Pyenv may require installing additional depencies.

If your distribution ships with a Python installation, do not attempt uninstall or update it. It may be an important part of your system.

For information on how to use Pyenv, see for example https://realpython.com/intro-to-pyenv/

The commands to install a virtual environment may look something like this:
```
cd server
python -m venv env
source env/Scripts/activate
```

----

## Getting started (macOS)

No instructions yet.

----

## Configuration file

Store environment-specific and sensitive settings in a config file. Keep the config file outside of any Git folders. Doing this minizes the risk of leaking credentials to Git or messing up other developers' configurations.

The config file uses normal Python syntax. See `default_settings.cfg` for more in-depth examples. Use uppercase keys for configuration options.

### Config file location

Set the path to the config file with an environment variable. Do not use quotes on Windows, even if the path contains spaces.

Windows:
```
set EPOTKU_SETTINGS=C:\path\to\epotku_settings.cfg
```

Linux / macOS:
```
export EPOTKU_SETTINGS=/path/to/epotku_settings.cfg
```

Protect the config file:
```
chmod 600 epotku_settings.cfg
```

### Secret key

To generate a secret key, run the following command:
```
python -c "import os; print(os.urandom(32).hex())"
```

Add the output of the previous command to `epotku_settings.cfg`:
```
SECRET_KEY = 'a85cdd9108e32dcd861b62aaa58f09a00810d2bbb635baeb19c7109d3344fbd4'
```

Generate and set another secret key for `JWT_SECRET_KEY`.

### More information

https://flask.palletsprojects.com/en/1.1.x/config/

----

# Other files

## Sample data

ePotku bundles a modified `Potku_v2.potku` as sample data for developing and testing. Logs, simulation files and miscellaneous other files have been deleted. The original file can be found [here](https://www.jyu.fi/science/en/physics/research/infrastructures/accelerator-laboratory/pelletron/potku/sample_data). The sample data is not licensed under GPL.

# Debris Removal Tool

@author: Yvan GARY, Gaëtan PIERRE, Côme OOSTERHOF, Myrtille MONCLIN, Tom SEMBLANET

## Introduction

Welcome to Debris Removal Tool's user guide ! On this page, you will find useful information in order to help you install the tool and use it. 

## Outline
 1. Description of Debris Removal Tool
 2. Installation procedure
 3. How to use Debris Removal Tool

# 1. Description of Debris Removal Tool

Debris Removal Tool (DRT) is a sample of python codes aiming to plan debris removal missions (consisting in multiple space rendezvous). For a given set of debris, DRT plans several missions, each one including 4 to 5 debris (this is a default setting that can be changed). 

It is focused on two major axes which are :
- the regroupement of the debris in order to plan optimal missions
- the plannification of the multiple rendezvous missions

## 2. Installation Procedure

In order to use DRT, you need to download the python codes in this repository. In order to do so, you need to have git installed on your machine.

### 2.1 Git installation

Git installation is different depending on you OS.

##### Linux - Debian / Mac OS

On Linux or Mac, open a command window and run the following command :

```bash
sudo apt install git
```

##### Windows

A git client application you can use: [TortoiseGit](https://tortoisegit.org/).
Another option: [SmartGit](https://www.syntevo.com/smartgit/).

### 2.2 Download DRT's content

Now that you have git installed on your machine, you are able to download DRT. First of all, open a command window and get to the location where you will install DRT. 

From here, run the following command on the terminal :

```bash
git clone https://github.com/TomSemblanet/PIE-COS05.git
```

This will create a folder name DRT containing python codes. Now, step into the folder :

```bash
cd DRT
```

The last step in now to download the dependencies needed to run DRT, which can be done by running the following command :

```bash
pip install -r requirements.txt
```

Debris Removal Tool in now ready to use !

## 3. How to use Debris Removal Tool

Make sure that you are into the DRT folder : 

```bash
cd DRT
```

From here, you can launch the program aiming to regroup the debris :

```bash
python -m regroupement.main
``` 

Note : The code is using by default a dataframe of 47 dangerous debris gathered by the IAF (International Astronautical Federation).




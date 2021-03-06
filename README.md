# OMCompiler

This repo is for local changes to the QSS solver and everything related to ultimately simulate the VHM in as short as possible time.

This is merely for the OMCompiler submodule.

## git remote -v

Ok, we need to pull the latest changes from OpenModelica's git repo (and sort out any merge conflicts that may arise manually) but push to this repo.


```
#!python

anton@qwerty:~/Distros/OpenModelica/OMCompiler$ git remote -v
origin	git@bitbucket.org:HealthQ/hq_openmodelica.git (fetch)
origin	git@bitbucket.org:HealthQ/hq_openmodelica.git (push)
upstream	https://openmodelica.org/git-readonly/OMCompiler.git (fetch)
upstream	https://openmodelica.org/git-readonly/OMCompiler.git (push)
```

To add another remote stream look at https://git-scm.com/book/no-nb/v1/Git-Branching-Remote-Branches

You will do something like


```
#!python

git remote add upstream2 git@bitbucket.org:HealthQ/hq_openmodelica.git
```




```
#!python

anton@qwerty:~/Distros/OpenModelica/.git/modules/OMC$ cat config
[core]
	repositoryformatversion = 0
	filemode = true
	bare = false
	logallrefupdates = true
	worktree = ../../../OMCompiler
[remote "origin"]
	url = git@bitbucket.org:HealthQ/hq_openmodelica.git
	fetch = +refs/heads/*:refs/remotes/origin/*
[branch "master"]
	remote = upstream
	merge = refs/heads/master
[submodule "3rdParty"]
	url = https://openmodelica.org/git-readonly/OMCompiler-3rdParty
[submodule "common"]
	url = https://openmodelica.org/git-readonly/OpenModelica-common
[remote "upstream"]
	url = https://openmodelica.org/git-readonly/OMCompiler.git
	fetch = +refs/heads/*:refs/remotes/upstream/*
```

## Operational Stuffs

1. You will also have to change the location of the wrtiting of the files which is found in the **LIQSS_write** function in the **liqss2Operations.c** file.

2. Also, the second derivatives in the **ddx** function must be compiled for the model under consideration.

3. The tolerances must be changed within the solver. The tolerances are named **tolerance** and **absTolerance**.

To build the model navigate to */OpenModelica/OMCompiler* and run *make* in the terminal.

# OMCompiler

This repo is for local changes to the QSS solver and everything related to ultimately simulate the VHM in as short as possible time.

This is merely for the OMCompiler submodule.

## git remote -v

Ok, we need to pull the latest changes from OpenModelica's git repo (and sort out any merge conflicts that may arise manually) but push to this repo.


```
#!python

anton@qwerty:~/Distros/OpenModelica/OMCompiler$ git remote -v
origin	https://antonpdv@bitbucket.org/antonpdv/omcompiler.git (fetch)
origin	https://antonpdv@bitbucket.org/antonpdv/omcompiler.git (push)
upstream	https://openmodelica.org/git-readonly/OMCompiler.git (fetch)
upstream	https://openmodelica.org/git-readonly/OMCompiler.git (push)
```

To add another remote stream look at https://git-scm.com/book/no-nb/v1/Git-Branching-Remote-Branches
 
You will do something like 


```
#!python

git remote add upstream2 https://openmodelica.org/git-readonly/OMCompiler.git
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
	url = https://antonpdv@bitbucket.org/antonpdv/omcompiler.git
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

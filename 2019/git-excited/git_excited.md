
class: center, middle

<img src="https://miro.medium.com/max/1400/1*mtsk3fQ_BRemFidhkel3dA.png" width="300" height="200" />

# Cupcakes and Coding: Git Excited for Git!

### Gabby Shvartsman

July 11th, 2019

image source: https://miro.medium.com/max/1400/1*mtsk3fQ_BRemFidhkel3dA.png

---

#### Image Citations:

All images I'll be using are from the [*Pro Git* online book](https://git-scm.com/book/en/v2) (unless cited otherwise)

---

## Start the Downloads

**Installing Git:**


* run `git` in your Terminal to see if you have it (if you don't, it'll prompt you to install it)


* if it just says it's not there and doesn't prompt installation, run `xcode-select --install`


**Setting Up Github:**


* make an account if you don't already have one (go to github.com and "Sign Up for Github")


* message Olga with your new or existing Github username, she'll add you to the czbiohub organization

---

## Things to Keep in Mind:

* it was originally designed by an advanced programmer for an advanced programmer


* so, using Git takes some time to get the hang of, don't worry!


* you'll probably mess up (good news: you can almost always recover it!!)


* best way to learn is to mess up and then fix it

---

## Useful Vocab:

* **repository:** a central location where data can be stored/managed
  * like a folder on your computer but now it is remote (outside of your own computer), other people can access
  * "repo"
  * Use: "Could you add the file to the Github repo?"
  
  
* **local vs. remote:**
  * local = on your computer, what you see with a file explorer
  * remote = on some other server (you don't have the files on your computer)
  * Use: "I have the file remotely but not locally"

---

class: center, middle

## What is Git?

<img src="https://cdn-images-1.medium.com/max/2600/1*y7D5jICmjzvxZP6z-5EtDg.png" width="300" height="200" />


image source: https://cdn-images-1.medium.com/max/2600/1*y7D5jICmjzvxZP6z-5EtDg.png

---

## Git is a fast & convenient *version control system*

### Version Control System:

* keeps track of versions of your files that you have saved in the past


* common example: **google docs "version history"**
  * (except doesn't automatically update the saved version, you have to do that yourself)
  * [docs version history example](https://docs.google.com/document/d/1UPuQWI4Xgx4oBQ0QL0BS5fqzf47hsW2QsXJi8wpuBdY/edit#)

---

## Types of Version Control Systems:

### Local:

* ex: you have different versions of your files in different folders on your computer so that you can go back to old versions
* super hard to maintain if you're working on a team

### Centralized:

* one central server holds files and you can only access certain files at a time to edit
* if server goes down, no one can work on it because you don't have all of the files saved locally

### *Distributed:*

* not just certain files, mirror of entire repository and all its history
* if server crashes, can still recover entire project

---

class: center, middle

## It's super common in software development

but also can use it for any file storage/saving!

---

<img src="https://slack-files2.s3-us-west-2.amazonaws.com/avatars/2017-12-19/288981919427_f45f04edd92902a96859_512.png" width="100" height="100" />

## Github:

* not the same as Git


* a display of all of your Git repos and all of the files in each repo


* can navigate like through file explorer


* useful features: 
  * past commits
  * git blame
  * branches
  * pull requests (PRs)

image source: https://slack-files2.s3-us-west-2.amazonaws.com/avatars/2017-12-19/288981919427_f45f04edd92902a96859_512.png

---

## Big ideas behind saving a file in Git:

![](https://git-scm.com/figures/18333fig0106-tn.png)

* all happening locally (no communication with a remote server yet)

---

## Common Local Git Commands

* `git add <filename>`: moves a file that you've been working on and want to save from your *working directory* into your *staging area*


* `git commit -m "<commit message>"`: moves everything from your staging area into your *local git repo* and takes down a message with it, so you can make a note to yourself about what you've changed in your file since the last time you committed any of your files


* at any point you can use `git status` to see which files are in which stages


* `git log` will let you see your commit history at any point

---

## Committing in Git

![](https://git-scm.com/figures/18333fig0105-tn.png)


* each commit you make is a "snapshot" of all the files you have in the directory

* files that have remained unchanged in the next commit will just have pointers to these older versions

* each commit can be referenced using a pointer

* the most recent commit you're using in your working directory has a pointer called HEAD

---

![](https://git-scm.com/figures/18333fig0302-tn.png)

---

## Basic Communication with a Remote Repo

* `git push`: allows you to take everything that you committed locally (any new changes you've saved) and update your *remote repo* with the changes from your *local repo*


* `git pull`: allows you to take any changes you've made remotely and update your *local repo* with the changes from your *remote repo*
  * made up of a `git fetch` (pulls down all info from the remote project) and a `git merge` (combines the two versions you're working with)
  * will try to merge your local and remote repos and will raise a `merge conflict` if something doesn't line up and you'll have to fix that before you can combine your local and remote repos
  * merging moves back to a common ancestor of the two versions you're trying to merge and creates a new snapshot with the changes from both (more on this when we get to branches)

---

## "Undoing" Things (you'll be doing a lot of this)

* `git reset HEAD` will let you reset your working directory to a past commit


* "unmodifying" a modified file: `git checkout --<file>` (if you realized you want the previous version of just that one file to be restored)


* google error messages! (if the git error messages don't tell you how to fix it)


* [good resource for getting out of a "git hole"](https://sethrobertson.github.io/GitFixUm/fixup.html)

---

class: middle

## Try out what we've covered so far before we move on with more commands

#### First time setup:

1) Make sure you've installed Git and you have your Github account set up
  * Open your Terminal and run `git --version` to check if you have Git

2) run `git config --global user.name "<your Git username>"` and `git config --global user.email <your Git email>` to save these so you don't need to pass them in every time you want to run a command

3) `git config --list` to make sure your configurations were updated

---

class: middle

#### Creating a remote repo + cloning it locally:

1) Go to the `repositories` tab in your Github profile

2) Select `New` and enter your `Repository Name`, check the `Initialize this repository with a README` box

3) Select `Clone or download`, select the `ssh` option and copy the link

4) In your Terminal, navigate (use cd commands) to whichever folder you want to import it into, enter `git clone <link>`

Note: you can also do `git clone <link> <dir-name>` to automatically make a new directory and clone your repo to it 

---

class: middle

#### Add + commit + push practice:

1) create a .txt file in your text editor in your cloned local repo

2) practice adding/committing/pushing it and seeing what it looks like in your remote repository once you've done that successfully

3) explore the commit history and git blame features

4) if you have time, try editing the file both remotely (through Github) and locally (in your text editor) before a `git pull` to cause a merge conflict and then follow the steps in the error messages to resolve it (edit same line remotely and locally)

---

## Branching

* default branch when you create a repo = master branch


* super useful for working on a project with multiple people where you're working on different parts of the code at the same time


* your commit history on your master branch stays put until you're ready to merge the commits from another branch


* don't want to both be working on master because you might run into merge conflicts, wouldn't want to experiment with a new feature that might break your existing code all on one single branch

---

## Real Life Analogy of the Usefulness of Branching

Imagine you're with your friends at Six Flags with your friends and there's a SUPER long line for this one ride you really want to go on. You've already been in line for 20 minutes but you're starting to get hungry and you also need sunscreen ASAP. You don't want to leave the line because that'll erase all the waiting you've done so far, so you send one friend to go get food and another friend to go get sunscreen while you decide to keep waiting in line.


You waiting in line are analogous to the default "master" branch that any repo starts on.

Each of your friends is a new "branch" in that repo.

The idea is that once your friends come back, you've made some progress in the line, so they join you at this new spot in line and also "merge" in their new changes (food and sunscreen in this case).

---

## How is a branch actually represented in Git?

![](https://git-scm.com/figures/18333fig0303-tn.png)


* it's a pointer to certain commit in your history (like the master branch I mentioned earlier)


* if you only have this one branch and you make a commit, the master pointer will move forward automatically as well


* think of this pointer as your family's position in line: as you make progress in line, your actual *position* wiwll automatically update with it to keep track


---

* Say you want to add a new "testing" feature to your project and you want to make this on a separate branch so you don't mess with the commit history in the default master branch


* `git checkout -b <branch-name>`: creates a new branch and "checks it out" or switches to it


* so you'd run `git checkout -b testing` and then your commit history looks like this:


![](https://git-scm.com/figures/18333fig0306-tn.png)


* git keeps track of which branch you currently have checked out using the "HEAD" pointer

---

* could do work on both of them and end up with something like this:

![](https://git-scm.com/figures/18333fig0309-tn.png)

---

could end up with something like this:

<img src="https://git-scm.com/figures/18333fig0316-tn.png" width="300" height="200" />


and when you merge you would get this:

<img src="https://git-scm.com/figures/18333fig0317-tn.png" width="300" height="200" />

---

## Pull Requests (PRs)

* could manually merge your master and testing branch, but Github has a more convenient way to do that: **pull requests**


* you can see all the changes that have been made in the merge, look at merge conflicts


* can set reviewers that have to look over and approve the request before merging (some repos require a reviewer before merging a PR, depends on how you set it up)

---

## Branching / PR Demo

1) Go to the cupcakes and coding repo: https://github.com/czbiohub/cupcakes


2) Clone the cupcakes repo: `git clone <link> dir-name`


3) Once you have it set up locally, make a new branch for yourself: `git checkout -b gfs/cupcakes-git` (replace initials with your initials/first name)
  * run `git branch` to look at all of your current branches and see your new one there

4) Navigate to 2019/git-excited/doggos directory through your command line


5) Find a dog picture you like on google and save it into that dir


6) Add and commit this new file to your local branch


7) Push your local branch to the remote repo


8) Submit a PR to merge it with master: set a description and also set the reviewers as Olga and the person next to you


9) Go through the process of approving the other person's PR

---

## Typical Project Workflow in Git

* lots of logical branches (not a new branch for every little thing but also don't end up doing most of the project in one single branch because then that defeats the purpose)


* pretty frequent commits


* merge branches pretty often to make sure the other people in your team are up to date with new features


* submitting PR's is important!! (especially in a project where you need approval before you merge)


* ["Github flow" runthrough](https://guides.github.com/introduction/flow/)

---

## Rebasing

* instead of having all of these branches in your history, you flatten out the commit history so it's very linear and straightforward


* easier to be clear about the timeline of your commits


* `git rebase <branch-name>`


* what it does: copies over all of your commits from the specified branch to your current branch and then adds your commits from your current branch on top of them

---

class: center

from this:

<img src="https://git-scm.com/figures/18333fig0331-tn.png" width="300" height="200" />

to this:

![](https://git-scm.com/figures/18333fig0335-tn.png)


---

## Additional Git Resources to Git Excited About

https://swcarpentry.github.io/git-novice/

http://think-like-a-git.net/

https://docs.google.com/document/d/1UPuQWI4Xgx4oBQ0QL0BS5fqzf47hsW2QsXJi8wpuBdY/edit#heading=h.6ti090mpqnma

https://git-scm.com/book/en/v2/Getting-Started-About-Version-Control  
(I personally used this one and really liked how in-depth it was! Was very good at showing the bigger picture of each command)

---

class: center, middle

### Fin :)

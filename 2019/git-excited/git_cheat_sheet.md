# Git Cheat-Sheet:

`git clone <link> <dir-name>`: cloning remote repo into new directory <dir-name>

`git init`: initialize a git repo in an existing local directory

`git status`: check what's up with your working directory/staging area

`git add <file-name>`: track and stage a file you modify

`git commit -m "<commit-message>"`: commit your edits to your local repo

`git commit --amend`: append any newly-staged files to the latest commit

`git push`: push new commits from your local repo to the remote repo

`<file-name>.gitignore`: if you don't want git to track something you're constantly changing

`git diff`: shows only things you've changed and haven't staged

`git diff --cached`: shows what you've changed including things that have already been staged

`git rm <file-name>`: remove a file

`git mv <file-from> <file-to>`: for moving and/or renaming files and tracking those changes

`git log`: view the commit history

`git last`: view last commit

'git push': push local branch commits to remote branch

`git branch <branch-name>`: make a new branch

`git checkout <branch-name>`: switch to an existing branch (moving HEAD to that branch)

`git checkout -b <branch-name>`: combines ^ those 2 steps

`git pull`: a `git fetch` and `git merge` from remote branch to local branch

`git rebase <branch-name>`: incorporate the commits from <branch-name> and append the commits from your current branch to the end of that

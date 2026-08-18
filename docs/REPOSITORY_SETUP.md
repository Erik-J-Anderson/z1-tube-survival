# Repository setup

These steps preserve the current working code before any directory reorganization.

## 1. Copy the scaffold files into the project root

Place these files beside the existing C++ source files:

```text
.gitignore
.gitattributes
CMakeLists.txt
README.md
```

## 2. Inspect what Git will track

From the project root:

```bash
git init
git add .
git status
```

Before committing, verify that large inputs and generated files such as these are **not**
staged:

```text
Z1+SP.dat
*.dump
*.data
*.log
*.csv
*.png
compiled executables
build/
```

If anything large is staged accidentally, do not delete the file from disk. Remove it only
from the Git index:

```bash
git restore --staged path/to/file
```

## 3. Make the baseline commit

```bash
git commit -m "Initial working tube survival implementation"
```

Optionally mark the known-working scientific baseline:

```bash
git tag -a v0.1-history-affine-baseline \
    -m "Working affine correction and history-dependent tube escape"
```

## 4. Create a private GitHub repository

Create an empty private repository on GitHub. Do not initialize it with another README,
`.gitignore`, or license because those files already exist locally.

Then connect the local repository to it. For SSH:

```bash
git branch -M main
git remote add origin git@github.com:YOUR_USERNAME/YOUR_REPOSITORY.git
git push -u origin main
git push origin v0.1-history-affine-baseline
```

For HTTPS, use the HTTPS repository URL instead:

```bash
git remote add origin https://github.com/YOUR_USERNAME/YOUR_REPOSITORY.git
```

## 5. Verify the remote

```bash
git remote -v
git status
git log --oneline --decorate -5
```

## Windows / WSL / cluster note

The included `.gitattributes` forces source files and shell scripts to LF line endings.
That is particularly useful when the same repository is edited on Windows/WSL and executed
on a Linux cluster.

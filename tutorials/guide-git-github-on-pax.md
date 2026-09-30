# Tutorial: Set Up Git and GitHub on Pax

Applied Bioinformatics - Tufts University Department of Biology

You will use Git and GitHub throughout the course to save work, share homework
repositories, and make your analyses easier to reproduce. This guide sets up GitHub
access from the Pax terminal once. After that, the everyday workflow is usually just:

```bash
git status
git add FILE
git commit -m "Describe the change"
git push
```

Work through the parts in order. Each part includes a check so you can tell whether the
setup worked before moving on.

## Before You Start

You need a GitHub account. If you do not have one, create it at <https://github.com>.
You may use whichever email you prefer; it does not have to be your Tufts address.

Open a terminal on Pax:

- OnDemand: <https://ondemand-prod.pax.tufts.edu/> > **Clusters > Tufts HPC Shell Access**
- SSH: `ssh your_utln@login-prod.pax.tufts.edu`

Everything in this guide is light and safe to run on a login node.

### One exception to the `$MYWORK` rule

In the labs, course files, data, and results should stay under your course work space,
such as `$COURSE/$USER`, rather than in your home directory. Git configuration and SSH
keys are the exception. They live in your home directory because that is where Git and
SSH look for them by design:

```text
$HOME/.gitconfig        your name, email, and Git preferences
$HOME/.ssh/             your SSH keys and SSH settings
```

These files are tiny settings files. The storage rule is about data and analysis output,
not login credentials.

## Part 1: Tell Git Who You Are

Git records a name and email on every commit. Set them once:

```bash
git config --global user.name "Your Full Name"
git config --global user.email "your.name@tufts.edu"
```

Two more settings reduce surprises later:

```bash
git config --global init.defaultBranch main
git config --global pull.rebase false
```

The first makes new repositories start on `main`, matching GitHub's default branch name.
The second tells Git how to combine work when a pull needs to merge changes.

Check your settings:

```bash
git config --global --list
```

You should see your name, your email, `init.defaultbranch=main`, and
`pull.rebase=false`.

> The email you use here appears in commits. If you would rather not publish a personal
> address, GitHub can give you a private `@users.noreply.github.com` address from
> **Settings > Emails**. Use that address in the `user.email` command instead.

## Part 2: Create an SSH Key

GitHub does not accept account passwords for Git operations. Use an SSH key so Pax can
prove to GitHub that it is allowed to push to your repositories.

An SSH key pair has two halves:

- The private key stays on Pax and should never be shared.
- The public key goes to GitHub.

Create a new key:

```bash
ssh-keygen -t ed25519 -C "your.name@tufts.edu"
```

You will be asked three things:

**Where to save it.** Press Enter to accept `~/.ssh/id_ed25519`.

**A passphrase.** For this course, pressing Enter twice for no passphrase is acceptable
because you already log in to Pax through Tufts authentication. If you prefer to set a
passphrase, Part 4 shows how to avoid retyping it repeatedly.

**Confirm the passphrase.** Press Enter again.

Check that both files exist:

```bash
ls -l ~/.ssh/id_ed25519*
```

You should see two files:

```text
~/.ssh/id_ed25519       private key; never share this
~/.ssh/id_ed25519.pub   public key; this is what you add to GitHub
```

### Permissions Matter

SSH refuses to use a private key if other people could read it. If you ever see a
message such as `Permissions 0644 for 'id_ed25519' are too open`, fix the permissions:

```bash
chmod 700 ~/.ssh
chmod 600 ~/.ssh/id_ed25519
chmod 644 ~/.ssh/id_ed25519.pub
```

> Never commit a private key. `id_ed25519` with no `.pub` ending is the secret half. If
> it ever ends up in a repository, ask for help and make a new key pair.

## Part 3: Add the Public Key to GitHub

Print the public key:

```bash
cat ~/.ssh/id_ed25519.pub
```

Copy the whole line. It starts with `ssh-ed25519` and ends with your email.

Then, on GitHub:

1. Click your avatar in the top-right corner and choose **Settings**
2. In the left sidebar, open **SSH and GPG keys**
3. Click **New SSH key**
4. Title: `Tufts Pax cluster`, or another name you will recognize
5. Key type: **Authentication Key**
6. Paste the public key line into **Key**
7. Click **Add SSH key**

If GitHub says the key is invalid, you probably copied only part of the line or included
an extra line break. Widen the terminal window and copy it again.

## Part 4: Test GitHub Login from Pax

Run:

```bash
ssh -T git@github.com
```

The first time, SSH may ask whether you trust GitHub's server:

```text
The authenticity of host 'github.com (...)' can't be established.
ED25519 key fingerprint is SHA256:+DiY3wvvV6TuJJhbpZisF/zLDA0zPMSvHdkr4UvCOqU.
Are you sure you want to continue connecting (yes/no/[fingerprint])?
```

That ED25519 fingerprint is GitHub's current public fingerprint. GitHub also lists its
fingerprints in the official documentation:
<https://docs.github.com/en/authentication/keeping-your-account-and-data-secure/githubs-ssh-key-fingerprints>

Type `yes`.

Success looks like this:

```text
Hi yourusername! You've successfully authenticated, but GitHub does not provide shell access.
```

The `does not provide shell access` part is not an error. It means GitHub accepted your
key for Git, but it will not open a normal shell session.

If the command hangs and then times out, go to Part 5. That is common on clusters and
does not mean you did anything wrong.

If it says `Permission denied (publickey)`, GitHub does not recognize your key. Re-check
that the public key in GitHub matches:

```bash
cat ~/.ssh/id_ed25519.pub
```

### If you set a passphrase

Start an SSH agent once per terminal session so you are not asked repeatedly:

```bash
eval "$(ssh-agent -s)"
ssh-add ~/.ssh/id_ed25519
```

## Part 5: Use Port 443 If Port 22 Is Blocked

Many clusters block outbound connections on port 22, the port SSH normally uses. GitHub
also accepts SSH connections on port 443, the HTTPS port, which is usually open.

Test the alternate route:

```bash
ssh -T -p 443 git@ssh.github.com
```

If you get the `Hi yourusername!` message, make that route permanent:

```bash
cat >> ~/.ssh/config <<'EOF'

Host github.com
    Hostname ssh.github.com
    Port 443
    User git
EOF

chmod 600 ~/.ssh/config
```

Now test the normal command again:

```bash
ssh -T git@github.com
```

It should succeed. From now on, standard GitHub SSH URLs such as
`git@github.com:YOURUSERNAME/appbio-hw.git` will quietly use port 443 from Pax.

## Part 6: Connect a Repository

### If you are creating a new homework repository

On GitHub, click **+** in the top-right corner, then **New repository**.

Recommended choices:

- Repository name: something clear, such as `appbio-hw`
- Visibility: **Private**, unless your instructor tells you otherwise
- Initialize with README, `.gitignore`, or license: leave these unchecked if you already have files on Pax

Then connect your Pax directory to that GitHub repository. Replace `YOURUSERNAME` and
`appbio-hw` as needed:

```bash
export COURSE=/cluster/tufts/bio_appbio
export MYWORK=$COURSE/$USER

cd "$MYWORK/appbio-hw"
git init
git branch -M main
git remote add origin git@github.com:YOURUSERNAME/appbio-hw.git
```

If the directory is already a Git repository, `git init` will not erase your files.

Before the first push, make at least one commit:

```bash
git status
git add README.md
git commit -m "Start homework repository"
git push -u origin main
```

If your first file is not named `README.md`, replace it with the file you want to commit.
The `-u origin main` part sets the default remote branch, so later pushes can be just
`git push`.

### If the repository already exists on GitHub

Use the SSH clone URL:

```bash
cd "$MYWORK"
git clone git@github.com:YOURUSERNAME/appbio-hw.git
cd appbio-hw
```

The SSH URL has this form:

```text
git@github.com:USERNAME/REPOSITORY.git
```

It uses a colon after `github.com` and does not start with `https://`.

## Part 7: Protect Your Repository Before You Push Data

It is much easier to keep large files and secrets out of Git than to remove them later.
Add a `.gitignore` before you start tracking analysis outputs:

```bash
cat > .gitignore <<'EOF'
# Sequencing data: never commit these
*.fastq
*.fastq.gz
*.fq.gz
*.bam
*.bai
*.sam
*.vcf.gz

# Generated directories
data/
results/
work/
.nextflow*

# Keys and secrets
id_ed25519
id_rsa
*.pem

# Editor and OS noise
.DS_Store
*~
EOF

git add .gitignore
git commit -m "Add gitignore for sequencing data and secrets"
```

Before you push, check size and status:

```bash
du -sh .
git status
```

Small repositories are expected. If the repository is hundreds of megabytes, pause and
check what is staged before pushing. GitHub warns above 50 MB per file and rejects files
above 100 MB.

## Everyday Git Cycle

Use this cycle whenever you work on a homework repository:

```bash
git status                       # see what changed
git add FILE                     # stage the file you want to keep
git commit -m "Describe change"  # save a snapshot
git push                         # send it to GitHub
```

Before turning in work, run:

```bash
git status
git log --oneline
du -sh .
```

Good signs:

- `git status` says the working tree is clean.
- `git log --oneline` shows a short history of meaningful commits.
- `du -sh .` reports a small repository, not large sequencing files.

## Troubleshooting

**`Permission denied (publickey)`** means GitHub does not recognize your key. Confirm
that the key exists with `ls ~/.ssh/id_ed25519`, confirm the public key is registered on
GitHub, and test again with `ssh -T git@github.com`.

**`Connection timed out`** usually means port 22 is blocked. Use Part 5.

**`Permissions 0644 for '...' are too open`** means SSH will not use the key until the
permissions are fixed. Use the `chmod` commands in Part 2.

**`Host key verification failed`** means the saved server identity does not match what
SSH sees now. Ask for help if you are unsure. If directed, remove the old entry and
reconnect:

```bash
ssh-keygen -R github.com
ssh -T git@github.com
```

**`Support for password authentication was removed`** usually means the remote URL uses
HTTPS instead of SSH. Switch the remote to SSH:

```bash
git remote -v
git remote set-url origin git@github.com:YOURUSERNAME/appbio-hw.git
git remote -v
```

**`Please tell me who you are`** means Part 1 was skipped or did not save correctly.
Re-run the `git config --global` commands.

**`fatal: not a git repository`** means you are not inside a Git repository. Use `cd` to
enter the repository directory, or run `git init` if this directory should become a
repository.

**`remote: error: File ... exceeds GitHub's file size limit`** means a large file was
committed. Ask for help rather than guessing; removing it correctly requires rewriting
Git history.

## Quick Reference

One-time setup:

```bash
git config --global user.name "Your Full Name"
git config --global user.email "your.name@tufts.edu"
git config --global init.defaultBranch main
git config --global pull.rebase false
ssh-keygen -t ed25519 -C "your.name@tufts.edu"
cat ~/.ssh/id_ed25519.pub
ssh -T git@github.com
```

If `ssh -T git@github.com` times out, add the port 443 block from Part 5 to
`~/.ssh/config`.

Everyday work:

```bash
git status
git add FILE
git commit -m "Describe change"
git push
```

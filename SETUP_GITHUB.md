# Setting Up TCRAssembler on GitHub

## Step 1: Create Repository on GitHub

1. Go to https://github.com/new
2. **Repository name**: `TCRAssembler`
3. **Description**: `A Python package for assembling T-cell receptor (TCR) sequences from IMGT templates or full contig sequences`
4. **Visibility**: Choose **Public**
5. **Important**: DO NOT check any boxes (no README, no .gitignore, no license - we already have these files)
6. Click **"Create repository"**

## Step 2: Initialize Git and Push (if not already a git repo)

After creating the repository on GitHub, run these commands in your terminal:

```bash
# Navigate to your project directory
cd /Users/yxu2/Documents/TCR_assembly

# Initialize git repository (if not already initialized)
git init

# Add all files (respecting .gitignore)
git add .

# Make initial commit
git commit -m "Initial commit: TCRAssembler - TCR Sequence Assembly Package"

# Add your GitHub remote (replace YOUR_USERNAME with your GitHub username)
git remote add origin https://github.com/YOUR_USERNAME/TCRAssembler.git

# Push to GitHub
git branch -M main
git push -u origin main
```

## Step 3: Verify

After pushing, verify your repository at:
`https://github.com/YOUR_USERNAME/TCRAssembler`

## Alternative: Using SSH (if you have SSH keys set up)

If you prefer using SSH instead of HTTPS:

```bash
git remote add origin git@github.com:YOUR_USERNAME/TCRAssembler.git
git push -u origin main
```

## Troubleshooting

### If you get "repository already exists" error:
- Check if you already initialized git: `git status`
- If yes, just add the remote and push:
  ```bash
  git remote add origin https://github.com/YOUR_USERNAME/TCRAssembler.git
  git push -u origin main
  ```

### If you need to change the remote URL:
```bash
git remote set-url origin https://github.com/YOUR_USERNAME/TCRAssembler.git
```

### If you get authentication errors:
- Make sure you're logged into GitHub
- You may need to use a Personal Access Token instead of password
- Or set up SSH keys for easier authentication


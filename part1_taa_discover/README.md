# Why is "Commit to main" Grayed Out in GitHub Desktop?

## Why is "Commit" Grayed Out?

- **No Changes Staged:**
  GitHub Desktop only allows committing if there are **staged changes** (files you’ve added, modified, or moved). If files are moved but not saved, or changes aren’t detected, the button will be grayed out.

- **No Files Selected:**
  You must **check the boxes** next to the files you want to commit.

- **Empty Commit Message:**
  GitHub Desktop requires a commit message, even if it’s brief.

- **Git Ignored Files:**
  Files listed in `.gitignore` won’t appear as changes.

---

## How to Fix It

### 1. Check for Unsaved Changes
- Save all files in your file explorer.
- In GitHub Desktop, go to the **"Repository"** menu and select **"Refresh"** (or press `F5`).

### 2. Select Files to Commit
- In the **"Changes"** tab, check the boxes next to the files or folders you want to commit.
- If changes don’t appear, ensure they were saved or moved correctly.

### 3. Add a Commit Message
- Enter a message in the **"Summary"** field (e.g., *"Move files to part2_antibody_design"*).

### 4. Check for Ignored Files
- Open `.gitignore` in your repository and confirm your files/folders aren’t listed.

### 5. Restart GitHub Desktop
- Closing and reopening the app can resolve temporary issues.

### 6. Check for Errors
- Look for error messages at the top of the GitHub Desktop window.

---

## Still Stuck?

### Manually Stage Files
Open a terminal in your repository and run:
```bash
git add .
git status

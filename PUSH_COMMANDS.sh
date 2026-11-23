#!/usr/bin/env bash
# ============================================================================
#  8 基因稿（P01）GitHub 上传命令 — 在 Git Bash 里逐段复制执行
#  仓库：gut-heart-8genes-AHG（新建，推荐）
#  包目录：F:/longxiawenjian/gaoxueyaqiguanzhou/Paper01_GutHeart_8genes/github_package_8genes
# ============================================================================

# ---------- 第 0 步：先确认 git 可用、身份已配置 ----------
git --version
git config --global user.name
git config --global user.email

# 若上面两条为空，先配一次（只需一次，之后永久生效）
git config --global user.name "Yan Wang"
git config --global user.email "18693323793@163.com"
git config --global core.autocrlf false
git config --global core.quotepath false


# ---------- 第 1 步：进入包目录 ----------
cd "F:/longxiawenjian/gaoxueyaqiguanzhou/Paper01_GutHeart_8genes/github_package_8genes"


# ---------- 第 2 步：初始化仓库并提交 ----------
git init -b main
git add -A
git status --short
git commit -m "Code and derived tables for the calibrated colocalization triage of eight gut-heart candidate genes"


# ---------- 第 3 步：关联远程仓库并推送 ----------
# 先在网页上建一个空仓库，名字 gut-heart-8genes-AHG，不要勾 README / .gitignore / licence
git remote add origin https://github.com/kongtonglangke/gut-heart-8genes-AHG.git
git push -u origin main


# ============================================================================
#  以下是备用与排错，平时不用跑
# ============================================================================

# 【改仓库名】如果网页上建的名字不是 gut-heart-8genes-AHG
# git remote set-url origin https://github.com/kongtonglangke/新仓库名.git

# 【已存在 origin】若提示 remote origin already exists
# git remote remove origin
# git remote add origin https://github.com/kongtonglangke/gut-heart-8genes-AHG.git

# 【认证失败】GitHub 不再接受账号密码，push 时：
#   用户名填 kongtonglangke
#   密码位置填 Personal Access Token（不是登录密码）
#   生成地址：https://github.com/settings/tokens  →  Generate new token (classic)
#   勾选 repo 权限，有效期选 90 天，生成后只显示一次，先粘到记事本
#   想免输入可用：git config --global credential.helper manager

# 【推到旧仓库】如果坚持用已有的 gut-heart-8genes-manuscript
# git remote add origin https://github.com/kongtonglangke/gut-heart-8genes-manuscript.git
# git pull --rebase origin main      # 旧仓库已有两次提交，先拉再推
# git push -u origin main

# 【想撤销提交但保留文件】
# git reset --soft HEAD~1

# 【看远程到底绑的哪个地址】
# git remote -v

# 【文件太大被拒】本包最大文件 0.1 MB，正常不会触发；若以后加大文件：
# git lfs install
# git lfs track "*.csv"


# ============================================================================
#  推完之后
# ============================================================================
# 1. 打开 https://github.com/kongtonglangke/gut-heart-8genes-AHG 确认文件齐全（约 75 个）
# 2. 点右上角 About 旁的齿轮，Description 填：
#    Code and derived tables for a three-layer genetic triage of eight gut-heart candidate genes in hypertension
# 3. Topics 填：colocalization  mendelian-randomization  eqtl  hypertension  gtex  finngen
# 4. 回到本机，把仓库地址写进论文的 Data Availability（若编辑要求 GitHub 链接）

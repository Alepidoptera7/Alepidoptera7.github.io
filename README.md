# Bioinformatics blog (Jekyll + GitHub Pages)

A personal site with a blog, projects, publications and an about page.
GitHub Pages builds it directly: no GitHub Actions, no build server, and only
plugins that GitHub Pages supports (`jekyll-feed`, `jekyll-seo-tag`,
`jekyll-sitemap`).

## Put it online (about 5 minutes)

1. On GitHub, create a **public** repository named exactly
   `yourusername.github.io` (with your real username).
2. Upload everything in this folder to the repository. Either:
   - drag the files into **Add file → Upload files** on the repo page, or
   - from a terminal:
     ```bash
     cd bioblog
     git init -b main
     git add .
     git commit -m "Initial site"
     git remote add origin https://github.com/Alepidoptera7/Alepidoptera7.github.io.git
     git push -u origin main
     ```
3. In the repo, go to **Settings → Pages**. Under *Build and deployment*, set
   **Source: Deploy from a branch**, **Branch: main**, folder **/ (root)**, and save.
4. After a minute or two the site is live at `https://yourusername.github.io`.
   The **Actions** tab shows each build and any errors.

**Using a different repo name** (e.g. `blog`)? The site then lives at
`https://yourusername.github.io/blog`, so set `baseurl: "/blog"` in `_config.yml`.

## Make it yours

| What                          | Where                         |
|-------------------------------|-------------------------------|
| Name, role, email, profiles   | `_config.yml`                 |
| Home page intro               | `index.md`                    |
| About page                    | `about.md`                    |
| Projects                      | `_data/projects.yml`          |
| Publications                  | `_data/publications.yml`      |
| Menu links                    | `_data/navigation.yml`        |
| Colours and fonts             | `assets/css/main.css` (top)   |

The three posts in `_posts/` and the entries in `_data/` are examples. Replace
or delete them. `_posts/2026-11-20-how-to-write-posts.md` shows every
formatting feature (code, equations, tables, figures, callouts).

## Write a post

Create `_posts/YYYY-MM-DD-short-title.md`:

```markdown
---
title: My first real post
description: One sentence for post lists and link previews.
tags: [rna-seq, tutorial]
math: false
---

Your text in Markdown...
```

Commit and push; the site rebuilds automatically. Posts dated in the future
are not published until that date.

## Preview on your computer (optional)

Needs Ruby 3.x and Bundler.

```bash
bundle install
bundle exec jekyll serve --livereload
# open http://localhost:4000
```

## Features

- Light and dark mode, following the reader's system setting
- Syntax highlighting for Python, R, Bash, YAML and more, with a Copy button on every code block
- LaTeX equations via MathJax, loaded only on posts with `math: true`
- Tags page, yearly archive, RSS feed (`/feed.xml`), sitemap and SEO/social preview tags
- Reading time and previous/next links on posts
- Publications list that underlines your name in author lists

## Layout

```
_config.yml          site settings
_data/               navigation, projects, publications
_includes/           head, header, footer, post list, helix graphic
_layouts/            default, home, page, post
_posts/              blog posts
assets/css/          main.css, syntax.css
assets/img/          images and favicon
index.md, about.md   home and about pages
blog.html, tags.html, projects.html, publications.html, 404.html
```

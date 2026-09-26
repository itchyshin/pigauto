#!/usr/bin/env python3
# Turns the rendered study.html into study_page.html for the claude.ai artifact viewer, which wraps a page in its own
# document skeleton: keep <title>, the fonts link and the <style> blocks at the top, then the body content.
import re
t = open("study.html").read()
styles = "\n".join(re.findall(r"<style[^>]*>.*?</style>", t, re.S))
body = re.search(r"<body[^>]*>(.*)</body>", t, re.S).group(1)
fonts = ('<link rel="stylesheet" href="https://fonts.googleapis.com/css2?family=IBM+Plex+Sans:wght@400;600'
         '&family=IBM+Plex+Mono&family=Source+Serif+4:ital,opsz,wght@0,8..60,400;0,8..60,600;1,8..60,400&display=swap">')
open("study_page.html", "w").write("<title>Rubin MI Study</title>\n" + fonts + "\n" + styles + "\n" + body)
print(len(styles), "chars of style;", len(body), "chars of body")

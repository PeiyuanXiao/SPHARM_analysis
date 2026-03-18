@echo off
echo Starting SPHARM Lab...
wsl -u root -e sh -c "mkdir -p /mnt/h && mount -t drvfs H: /mnt/h"
docker start spharm-processor
start http://localhost:8787
echo Done! RStudio is opening in your browser.
pause
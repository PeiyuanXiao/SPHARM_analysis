@echo off
echo Starting SPHARM Lab...

:: 检查 Docker 是否在运行
docker info >nul 2>&1
if errorlevel 1 (
    echo [ERROR] Docker is not running. Please start Docker Desktop first.
    pause
    exit /b 1
)

:: 挂载 H 盘（已挂载则跳过）
wsl -u root -e sh -c "mountpoint -q /mnt/h || (mkdir -p /mnt/h && mount -t drvfs H: /mnt/h)"

:: 用 docker compose 启动容器
docker compose -f H:\SPHARM_analysis\docker-compose.yml up -d
if errorlevel 1 (
    echo [ERROR] Failed to start spharm-processor.
    pause
    exit /b 1
)

:: 等待 RStudio 就绪
echo Waiting for RStudio to be ready...
:wait_loop
curl -s http://localhost:8787 >nul 2>&1
if errorlevel 1 (
    timeout /t 2 /nobreak >nul
    goto wait_loop
)

start http://localhost:8787
echo Done! RStudio is opening in your browser.
pause
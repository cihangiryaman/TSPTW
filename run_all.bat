@echo off
echo ============================================
echo   TSPTW Solver - Build and Run
echo ============================================
echo.

echo [1/4] Compiling...
gcc -O3 -march=native -Wall -o tsptw.exe tsptw.c -lm
if errorlevel 1 (
    echo COMPILATION FAILED
    pause
    exit /b 1
)
echo Compilation successful.
echo.

echo [2/4] Running example-input-1.txt (76 cities)...
tsptw.exe example-input-1.txt output-1.txt 30
echo.

echo [3/4] Running example-input-2.txt (280 cities)...
tsptw.exe example-input-2.txt output-2.txt 30
echo.

echo [4/4] Running example-input-3.txt (~15000 cities)...
tsptw.exe example-input-3.txt output-3.txt 60
echo.

echo ============================================
echo   All runs complete!
echo ============================================
pause

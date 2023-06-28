rem     Uses cl.exe to create plain-knit.exe
rem     Ensure that cl.exe is on your PATH 
rem     (e.g., by using a Visual Studio Developer Command Prompt)

cl.exe /O2 /Wall /Fe"plain-knit.exe" plain-knit.c
del plain-knit.obj

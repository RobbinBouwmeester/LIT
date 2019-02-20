Set oShell = CreateObject ("Wscript.Shell") 
Dim strArgs
strArgs = "cmd /c LIT.bat"
oShell.Run strArgs, 0, false
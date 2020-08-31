# system(paste("/usr/lib/x86_64-linux-gnu/modulecmd.tcl python","avail"),intern=T)
# 
system(paste("/usr/lib/x86_64-linux-gnu/modulecmd.tcl python","load bioinfo-tools InfoMap"),intern=T)

module <- function(Arguments, moduleCmd="/usr/lib/x86_64-linux-gnu/modulecmd.tcl") { 
  pythonCmds <- system(paste(moduleCmd,"python",Arguments),intern=T) 
  print(pythonCmds)
  #eval(cmds) 
    # Check if all python commands are recognizable
  validPythonCmd <- grepl("os\\.chdir\\('([^']*)'\\)",pythonCmds) |
    grepl("os\\.environ\\['([^']*)'] = '([^']*)'",pythonCmds) |
    grepl("del os\\.environ\\['([^']*)'\\]",pythonCmds)
  if( !all(validPythonCmd) ){
    stop("modulecmd returned unknown command(s):\n", paste(pythonCmds[!validPythonCmd],collapse = "\n"))
  }
  
  # convert python commands to R commands
  RCmds <- sub("os\\.chdir\\('([^']*)'\\)","setwd(dir = '\\1')",pythonCmds,perl=T)
  RCmds <- sub("os\\.environ\\['([^']*)'] = '([^']*)'","Sys.setenv('\\1' = '\\2')",RCmds,perl=T)
  RCmds <- sub("del os\\.environ\\['([^']*)'\\]","Sys.unsetenv('\\1')",RCmds,perl=T)
  
  # execute R commands
  print(RCmds)
  invisible( eval( parse(text = RCmds) ) )
}
module("avail")
module("load bioinfo-tools InfoMap")

system(paste("/usr/bin/tclsh",
             "/usr/lib/x86_64-linux-gnu/modulecmd.tcl", 
             #"python","list"),intern=T)
              "python","load", "bioinfo-tools InfoMap"),intern=T)

bla <- system(paste("/usr/lib/x86_64-linux-gnu/modulecmd.tcl", 
             "python","avail"),intern=T)

pythonCmds <- system(paste("/usr/lib/x86_64-linux-gnu/modulecmd.tcl", 
             "python","load", "bioinfo-tools InfoMap"),intern=T)
validPythonCmd <- grepl("os\\.chdir\\('([^']*)'\\)",pythonCmds) |
  grepl("os\\.environ\\['([^']*)'] = '([^']*)'",pythonCmds) |
  grepl("del os\\.environ\\['([^']*)'\\]",pythonCmds)

if( !all(validPythonCmd) ){
  stop("modulecmd returned unknown command(s):\n", paste(pythonCmds[!validPythonCmd],collapse = "\n"))
}

RCmds <- sub("os\\.chdir\\('([^']*)'\\)","setwd(dir = '\\1')",pythonCmds[2:7],perl=T)
RCmds <- sub("os\\.environ\\['([^']*)'] = '([^']*)'","Sys.setenv('\\1' = '\\2')",RCmds,perl=T)
RCmds <- sub("del os\\.environ\\['([^']*)'\\]","Sys.unsetenv('\\1')",RCmds,perl=T)

invisible( eval( parse(text = RCmds) ) )

system(paste("Infomap"),intern=T)

pythonCmds <- system(paste("/usr/lib/x86_64-linux-gnu/modulecmd.tcl", 
                           "r","load", "bioinfo-tools FastQC"),intern=T)
( eval( parse(text = pythonCmds) ) )

module <- function(Arguments, moduleCmd="/usr/lib/x86_64-linux-gnu/modulecmd.tcl") { 
  rCmds <- system(paste(moduleCmd, "r", Arguments),intern=T) 
  ( eval( parse(text = rCmds) ) )
}


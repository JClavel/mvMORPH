## Echo a message when loading the package
 
msg= c("##",paste0(" On est deja grands : .................")," ##\n",
       "##",rep(" ",14),paste0(" ,--.  ,--."),rep(" ",14)," ##\n",
       "##",rep(" ",14),paste0("/   | /    ","\\"),rep(" ",14),"##\n",
       "##",rep(" ",14),paste0("\`","|  ||  ()  |"),rep(" ",12)," ##\n",
       "##",rep(" ",14),paste0(" |  | ","\\","    /"),rep(" ",14),"##\n",
       "##",rep(" ",14),paste0(" \`","--","\'","  ","\`","--","\'"),rep(" ",14)," ##\n",
       "##",rep(' ',16),'_','|','____','|','_',rep(' ',15),' ##\n',
       "##",rep(' ',15),'|',rep(' ',8),'|',rep(' ',14),' ##\n',
       "##",rep(' ',12),'+',rep("-",14),'+',rep(' ',11),' ##\n',
       "##",rep(' ',12),'|',rep(' ',14),'|',rep(' ',11),' ##\n',
       "##",rep(' ',8),'+',rep("-",22),'+',rep(' ',7),' ##\n',
      # "##",rep(' ',39),' ##\n',
       "##",rep(' ',8),'|',rep(' ',22),'|',rep(' ',7),' ##\n',
       
       #"##",paste0("        |                    |        "),"  ##\n",
       "##",paste0("   + ------------------------------ + "),"  ##\n",
       "##",paste0("            __  __  ___  ___ ___ _  _ "),"  ##\n",
       "##",paste0("   _ ____ _|  ","\\","/  |/ _ ","\\","| _ ","\\"," _ ","\\ || |"),"  ##\n",
       "##",paste0("  | '  ","\\"," V / |","\\","/| | (_) |   /  _/ __ |"),"  ##\n",
       "##",paste0("  |_|_|_","\\","_/|_|  |_|","\\___/|_|_","\\_| |_||_|"),"  ##\n",
       "##",paste0("... et ce n'est que le debut :)       "),"  ##\n")

.onAttach <- function(...) {
    # echo output to screen
    packageStartupMessage("##\n## mvMORPH package (1.2.2)")
    packageStartupMessage(msg)
    packageStartupMessage("## Multivariate evolutionary models")
    packageStartupMessage("##\n## See the tutorials: browseVignettes(\"mvMORPH\")")
    packageStartupMessage("##\n## To cite package 'mvMORPH': citation(\"mvMORPH\")\n##")

}

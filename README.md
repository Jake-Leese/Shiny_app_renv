Shiny app repository for exploring Chick Neural Plate Border scATAC-seq data. Also includes renv for reproducing the compatible R package versions.

Global.R, ui.R and server.R contain the original script set up to work with ShinyServer

Dynamic_Shiny_Script.R provides a workaround for launching the app with minimal compute resources (< 8GB mem) by dynamically loading and unloading each ArchR object as and when it is called by the application

To set up Shiny Server:
- After installing Shiny Server, clone the repo into "/srv/shiny-server/"
- Open directory in R and use renv::restore(), followed by renv::isolate() to load R packages and remove reliance on the renv cache
- Remove Dynamic_Shiny_Script.R or use it to replace contents of Global.R, ui.R and server.R, depending on resources available
- Ensure rwx permissions are granted to Shiny Server and then the app should be ready to use

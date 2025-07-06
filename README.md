renv repository for Shiny app.

global.R, server.R, ui.R and custom_functions.R are sufficient for shiny server to run the application
Dynamic_Shiny_Script.R and Dynamic_Shiny_optimised.R offer less resource intense options for running the app, in that they dynamically unload and load ArchR objects as they are called by the shiny dashboard.

To set up Shiny server:
- Ensure R is installed along with Shiny and renv
- Install Shiny Server and pull Shiny_app_renv into /srv/shiny-server/
- Launch R from this directory and submit following commands:
- renv::restore()
- renv::isolate()
- Depending on computational resources, delete Dynamic_Shiny* scripts, or use them to replace contents of global.R, ui.R and server.R
- The application should be available to use at your_ip:3838/Shiny_app_renv/
- See /var/log/shiny-server/ to view log files and troubleshoot any issues with the application

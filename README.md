# FUGU-MS
Filtering utility for grouping untargeted mass spectrometry datasets.

![FUGU-MS Logo](Application/www/assets/Logo.png)

## FUGU-App
This application is available on our server: https://lewisresearchgroup.shinyapps.io/fugu-ms/

<!-- > [FUGU-MS Web Application](https://lewisresearchgroup.shinyapps.io/fugu-ms/).  -->

![FUGU-App Landing Page](Application/www/assets/FUGU-AppLandingPage.png)

To deploy the application locally, the user must have R and the following packages installed.

```
install.packages("shiny")
install.packages("shinyWidgets")
install.packages("shinybusy")
install.packages("plotly")
install.packages("ggplot2")
install.packages("tidyr")
install.packages("dplyr")
```

To run the app, simply open an R console, set the directory to the project folder, and deploy the app.

```
setwd("~/FUGO-MS/Application")
runApp()
```

This will launch a local instance of the app on the user's machine.

# FUGU-MS
Filtering utility for grouping untargeted mass spectrometry datasets.

![FUGU-MS Logo](Application/assets/Logo.png)

## FUGU-App
This application is available on our server: https://lewisresearchgroup.shinyapps.io/fugo-ms/

<!-- > [FUGU-MS Web Application](https://lewisresearchgroup.shinyapps.io/fugo-ms/).  -->

![FUGU-App Landing Page](Application/assets/FUGU-AppLandingPage.png)

To deploy the application locally, the user must have R and the following packages installed.

```
install.packages('shiny')
install.packages("shinyWidgets")
install.packages("shinybusy")
install.packages('plotly')
```

To run the app, simply open an R console, set the directory to the project folder, and deploy the app.

```
setwd("~/FUGO-MS")
runApp()
```

This will launch a local instance of the app on the user's machine.

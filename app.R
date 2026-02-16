# --- REQUIRED LIBRARIES ---
if (!require("pacman")) install.packages("pacman")
pacman::p_load(shiny, shinydashboard, ggplot2, randomForest, 
               xgboost, e1071, rpart, rpart.plot, writexl, dplyr, tidyr)

options(shiny.maxRequestSize = 100 * 1024^2)

# --- UI SECTION ---
ui <- dashboardPage(
  skin = "blue",
  dashboardHeader(title = "InSe Unified Lab"),
  dashboardSidebar(
    sidebarMenu(
      menuItem("1. Data Loading", tabName = "upload", icon = icon("upload")),
      menuItem("2. Unified Metrics", tabName = "metrics", icon = icon("table")),
      menuItem("3. Unified Importance", tabName = "importance", icon = icon("balance-scale")),
      menuItem("4. Comparative Parity", tabName = "parity", icon = icon("chart-line")),
      menuItem("5. Error Analysis", tabName = "errors", icon = icon("bug"))
    )
  ),
  dashboardBody(
    tabItems(
      tabItem(tabName = "upload",
              fluidRow(
                box(title = "Nanosonde (ns) Data", width = 6, status = "primary", solidHeader = TRUE,
                    fileInput("files_ns", "Upload ns Files", multiple = TRUE)),
                box(title = "Picosonde (ps) Data", width = 6, status = "info", solidHeader = TRUE,
                    fileInput("files_ps", "Upload ps Files", multiple = TRUE))
              ),
              fluidRow(
                box(title = "Processing Configurations", width = 12,
                    column(3, numericInput("skip", "Skip Rows:", 0)),
                    column(3, selectInput("sep", "Separator:", choices = c("Tab" = "\t", "Space" = "", "Comma" = ","))),
                    column(3, numericInput("seed", "Seed:", 123)),
                    column(3, sliderInput("train_split", "Train Ratio %:", 50, 95, 80)))
              )),
      
      tabItem(tabName = "metrics",
              fluidRow(
                box(title = "Combined Performance (All Models)", width = 12, 
                    tableOutput("unified_table"), downloadButton("dl_excel", "Export Data (Excel)")),
                box(title = "R-Squared Benchmark", width = 12, 
                    plotOutput("unified_metrics_plot"), downloadButton("dl_metrics_plot", "Download (600 DPI)"))
              )),
      
      tabItem(tabName = "importance",
              fluidRow(
                box(title = "Feature Importance: ns vs ps", width = 12,
                    plotOutput("unified_imp_plot", height = "500px"),
                    downloadButton("dl_imp_plot", "Download (600 DPI)"))
              )),
      
      tabItem(tabName = "parity",
              fluidRow(
                box(title = "Parity Plot: ns", width = 6, plotOutput("parity_ns"), downloadButton("dl_p_ns", "Download ns")),
                box(title = "Parity Plot: ps", width = 6, plotOutput("parity_ps"), downloadButton("dl_p_ps", "Download ps"))
              )),
      
      tabItem(tabName = "errors",
              fluidRow(
                box(title = "Residual Density", width = 12,
                    plotOutput("unified_resid_plot"),
                    downloadButton("dl_resid_plot", "Download (600 DPI)"))
              ))
    )
  )
)

# --- SERVER SECTION ---
server <- function(input, output, session) {
  pub_theme <- theme_bw(base_size = 15) + theme(legend.position = "top")
  # n₂ (nonlinear refractive index) hesaplama fonksiyonu
  calculate_n2 <- function(delta_T_pv, S, I0, lambda = 1064e-9, n0 = 2.5, L = 1e-3, alpha0 = 0) {
    # delta_T_pv: Peak-valley transmittance difference
    # S: Aperture linear transmittance (S = 1 - exp(-2*ra^2/wa^2))
    # I0: On-axis irradiance at focus (W/m^2)
    # lambda: Wavelength (m)
    # n0: Linear refractive index (typical for InSe ~2.5)
    # L: Sample thickness (m)
    # alpha0: Linear absorption coefficient (m^-1)
    
    # Effective length
    Leff <- if(alpha0 == 0) L else (1 - exp(-alpha0 * L)) / alpha0
    
    # On-axis phase shift
    delta_phi0 <- delta_T_pv / (0.406 * (1 - S)^0.25)
    
    # n2 in m^2/W
    k <- 2 * pi / lambda
    n2_m2W <- delta_phi0 / (k * I0 * Leff)
    
    # Convert to esu (1 m^2/W = (c*n0)/(40*pi) esu)
    c <- 3e8  # speed of light in m/s
    n2_esu <- n2_m2W * (c * n0) / (40 * pi)
    
    return(list(
      delta_phi0 = delta_phi0,
      n2_m2W = n2_m2W,
      n2_esu = n2_esu,
      Leff = Leff
    ))
  }
  
  # β (nonlinear absorption coefficient) hesaplama fonksiyonu
  calculate_beta <- function(T_z, I0, Leff) {
    # T_z: Normalized transmittance at focus (open aperture)
    # I0: On-axis irradiance at focus (W/m^2)
    # Leff: Effective length (m)
    
    # β in m/W
    beta_mW <- (2 * sqrt(2) * (1 - T_z)) / (I0 * Leff)
    
    # Convert to cm/GW (1 m/W = 10^9 cm/GW)
    beta_cmGW <- beta_mW * 1e9 * 100  # m/W -> cm/GW
    
    return(list(
      beta_mW = beta_mW,
      beta_cmGW = beta_cmGW
    ))
  }
  
  # I0 hesaplama (odakta ışın şiddeti)
  calculate_I0 <- function(P_in, f = 0.2, beam_diameter = 0.005, lambda = 1064e-9) {
    # P_in: Input power (W)
    # f: Focal length (m)
    # beam_diameter: Beam diameter before lens (m)
    # lambda: Wavelength (m)
    
    # Beam waist at focus (Gaussian beam)
    w0 <- (4 * lambda * f) / (pi * beam_diameter)
    
    # Beam area at focus
    A <- pi * w0^2 / 2  # for Gaussian beam
    
    # Intensity (W/m^2)
    I0 <- (2 * P_in) / A  # peak intensity for Gaussian
    
    return(list(
      w0 = w0,
      I0 = I0
    ))
  }
  # --- Robust Data Processing with Encoding Fix ---
  process_raw <- function(files) {
    req(files)
    file_list <- list()
    for(i in 1:nrow(files)) {
      temp <- read.table(files$datapath[i], skip = input$skip, sep = input$sep, fill = TRUE, stringsAsFactors = FALSE)
      # Force columns
      W <- as.numeric(as.character(temp[,1]))
      I <- as.numeric(as.character(temp[,2]))
      
      # Safer Encoding Handling
      fname <- files$name[i]
      fname_clean <- iconv(fname, from = "", to = "UTF-8", sub = "byte")
      # Extract numeric value safely
      d_val_str <- gsub("[^0-9,.]", "", fname_clean)
      d_val <- as.numeric(gsub(",", ".", d_val_str))
      if(is.na(d_val)) d_val <- 0
      
      clean_df <- data.frame(Wavelength = W, Intensity = I, Doping = d_val)
      file_list[[i]] <- clean_df[complete.cases(clean_df), ]
    }
    full_df <- do.call(rbind, file_list)
    set.seed(input$seed)
    idx <- sample(1:nrow(full_df), floor(nrow(full_df) * (input$train_split/100)))
    list(train = full_df[idx, ], test = full_df[-idx, ])
  }
  
  data_ns <- reactive({ process_raw(input$files_ns) })
  data_ps <- reactive({ process_raw(input$files_ps) })
  
  # --- Unified Model Engine ---
  engine <- reactive({
    req(data_ns(), data_ps())
    
    r2_f <- function(a, p) { 1 - (sum((a-p)^2)/sum((a-mean(a))^2)) }
    rmse_f <- function(a, p) { sqrt(mean((a - p)^2)) }
    mae_f <- function(a, p) { mean(abs(a - p)) }
    
    get_all <- function(ds, tag) {
      set.seed(input$seed)
      m_lr <- lm(Intensity ~ Wavelength + Doping, data = ds$train)
      m_dt <- rpart(Intensity ~ Wavelength + Doping, data = ds$train)
      m_rf <- randomForest(Intensity ~ Wavelength + Doping, data = ds$train)
      m_svm <- svm(Intensity ~ Wavelength + Doping, data = ds$train)
      m_xgb <- xgboost(data = as.matrix(ds$train[,c("Wavelength","Doping")]), label = ds$train$Intensity, nrounds=50, verbose=0)
      
      mods <- list(LR=m_lr, DT=m_dt, RF=m_rf, SVM=m_svm, XGB=m_xgb)
      res_table <- data.frame()
      for(n in names(mods)) {
        p <- if(n=="XGB") predict(mods[[n]], as.matrix(ds$test[,c("Wavelength","Doping")])) else predict(mods[[n]], ds$test)
        act <- ds$test$Intensity
        res_table <- rbind(res_table, data.frame(Model=n, TimeScale=tag, R2=r2_f(act, p), RMSE=rmse_f(act, p), MAE=mae_f(act, p)))
      }
      list(table = res_table, rf = m_rf)
    }
    
    list(ns = get_all(data_ns(), "ns"), ps = get_all(data_ps(), "ps"))
  })
  
  # --- Plots as Reactives (For Stability) ---
  plot_r2 <- reactive({
    df <- rbind(engine()$ns$table, engine()$ps$table)
    ggplot(df, aes(Model, R2, fill=TimeScale)) + geom_bar(stat="identity", position="dodge") + 
      scale_fill_manual(values=c("ns"="#E69F00", "ps"="#56B4E9")) + pub_theme
  })
  
  plot_imp <- reactive({
    i_ns <- as.data.frame(importance(engine()$ns$rf)); i_ps <- as.data.frame(importance(engine()$ps$rf))
    i_ns$TimeScale <- "ns"; i_ps$TimeScale <- "ps"; i_ns$Feature <- rownames(i_ns); i_ps$Feature <- rownames(i_ps)
    ggplot(rbind(i_ns, i_ps), aes(Feature, IncNodePurity, fill=TimeScale)) + 
      geom_bar(stat="identity", position="dodge") + scale_fill_manual(values=c("ns"="#E69F00", "ps"="#56B4E9")) + pub_theme
  })
  
  plot_resid <- reactive({
    r_ns <- data.frame(err = data_ns()$test$Intensity - predict(engine()$ns$rf, data_ns()$test), ts = "ns")
    r_ps <- data.frame(err = data_ps()$test$Intensity - predict(engine()$ps$rf, data_ps()$test), ts = "ps")
    ggplot(rbind(r_ns, r_ps), aes(err, fill=ts)) + geom_density(alpha=0.5) + 
      scale_fill_manual(values=c("ns"="#E69F00", "ps"="#56B4E9")) + pub_theme + labs(x="Residual")
  })
  
  # --- RENDERING & EXPORTS ---
  output$unified_table <- renderTable({ rbind(engine()$ns$table, engine()$ps$table) })
  output$unified_metrics_plot <- renderPlot({ plot_r2() })
  output$unified_imp_plot <- renderPlot({ plot_imp() })
  output$unified_resid_plot <- renderPlot({ plot_resid() })
  
  output$parity_ns <- renderPlot({ 
    ds <- data_ns(); p <- predict(engine()$ns$rf, ds$test)
    ggplot(data.frame(a=ds$test$Intensity, p), aes(a,p)) + geom_point(color="#E69F00", alpha=0.4) + geom_abline() + pub_theme + labs(title="ns Scale")
  })
  output$parity_ps <- renderPlot({ 
    ds <- data_ps(); p <- predict(engine()$ps$rf, ds$test)
    ggplot(data.frame(a=ds$test$Intensity, p), aes(a,p)) + geom_point(color="#56B4E9", alpha=0.4) + geom_abline() + pub_theme + labs(title="ps Scale")
  })
  
  # --- DOWNLOAD HANDLERS ---
  output$dl_excel <- downloadHandler(filename="InSe_Unified_Metrics.xlsx", content=function(f) write_xlsx(rbind(engine()$ns$table, engine()$ps$table), f))
  output$dl_metrics_plot <- downloadHandler(filename="R2_Comparison_600DPI.png", content=function(f) ggsave(f, plot_r2(), dpi=600))
  output$dl_imp_plot <- downloadHandler(filename="Importance_Comparison_600DPI.png", content=function(f) ggsave(f, plot_imp(), dpi=600))
  output$dl_resid_plot <- downloadHandler(filename="Residuals_Comparison_600DPI.png", content=function(f) ggsave(f, plot_resid(), dpi=600))
  output$dl_p_ns <- downloadHandler(filename="Parity_ns_600DPI.png", content=function(f) {
    ds <- data_ns(); p <- predict(engine()$ns$rf, ds$test)
    ggsave(f, ggplot(data.frame(a=ds$test$Intensity, p), aes(a,p)) + geom_point(color="#E69F00", alpha=0.4) + geom_abline() + pub_theme, dpi=600)
  })
  output$dl_p_ps <- downloadHandler(filename="Parity_ps_600DPI.png", content=function(f) {
    ds <- data_ps(); p <- predict(engine()$ps$rf, ds$test)
    ggsave(f, ggplot(data.frame(a=ds$test$Intensity, p), aes(a,p)) + geom_point(color="#56B4E9", alpha=0.4) + geom_abline() + pub_theme, dpi=600)
  })
}

shinyApp(ui, server)

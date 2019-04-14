#' Display the graphical user interface.
#' @useDynLib AbdominalEIT
#' @export
eit <- function() {
    # Sys.setenv('http_proxy'='http://proxy.sdd.kao.co.jp:8080')
    options(device = "x11")
    # options(device='quartz') imaging of EIT data Start point
    rm(list=ls())
    prep <- function(obj) {
        # browser()
        aniso <- gtkComboBoxGetActive(sw1)
        preProcess(anisotropy = aniso)
    }
    recon <- function(obj) {
        # browser() library(rgl)
        aniso <- gtkComboBoxGetActive(sw1)
        resEnh <- gtkComboBoxGetActive(sw4) + 1
        pattern <- gtkComboBoxGetActive(sw5)
        imageReconstruction(anisotropy = aniso, resEnh = resEnh, pattern = pattern)
    }
    sim <- function(obj) {
        # browser()
        aniso <- gtkComboBoxGetActive(sw1)
        model <- gtkComboBoxGetActive(sw2) + 1
        interval <- gtkComboBoxGetActive(sw3) + 1
        resEnh <- gtkComboBoxGetActive(sw4) + 1
        pattern <- gtkComboBoxGetActive(sw5)
        simulationStudy(model = model, interval = interval, resEnh = resEnh, anisotropy = aniso, pattern = pattern)
    }
    wn <- length(gtkWindowListToplevels())
    wn
    if (wn == 0) {
        win <- gtkWindowNew(show = FALSE)
        # gtkWindowSetResizable(win,FALSE)
        win["default-width"] <- 385
        win["default-height"] <- 290
        win["title"] <- "Main Menu"
        button1 <- gtkButtonNewWithLabel("Pre processing")
        button2 <- gtkButtonNewWithLabel("Image reconstruction")
        button3 <- gtkButtonNewWithLabel("Post processing")
        button4 <- gtkButtonNewWithLabel("Simulation study")
        sw1 <- gtkComboBoxNewText(show = TRUE)
        gtkComboBoxAppendText(sw1, "without considering anisotropy")
        gtkComboBoxAppendText(sw1, "with considering anisotropy")
        sw2 <- gtkComboBoxNewText(show = TRUE)
        gtkComboBoxAppendText(sw2, "Model 1: Notch, 0.85; otherwise, 0.75")
        gtkComboBoxAppendText(sw2, "Model 2: Notch, 0.75; otherwise, 0.85")
        gtkComboBoxAppendText(sw2, "Model 3: Notch, 0.90; otherwise, 0.70")
        gtkComboBoxAppendText(sw2, "Model 4: Notch, 0.70; otherwise, 0.90")
        gtkComboBoxAppendText(sw2, "Model 5: Notch, 0.85; otherwise, 0.75; abdominal muscle")
        gtkComboBoxAppendText(sw2, "Model 6: Notch, 0.75; otherwise, 0.85; abdominal muscle")
        gtkComboBoxAppendText(sw2, "Model 7: Notch, 0.90; otherwise, 0.70; abdominal muscle")
        gtkComboBoxAppendText(sw2, "Model 8: Notch, 0.70; otherwise, 0.90; abdominal muscle")
        gtkComboBoxAppendText(sw2, "Model 9: Notch, 0.90; otherwise, 0.70; muscle w/ aniso")
        gtkComboBoxAppendText(sw2, "Model 10: Notch, 0.70; otherwise, 0.90; muscle w/ aniso")
        gtkComboBoxAppendText(sw2, "Model 11: Realistic model Type 1 - Small")
        gtkComboBoxAppendText(sw2, "Model 12: Realistic model Type 1 - Medium")
        gtkComboBoxAppendText(sw2, "Model 13: Realistic model Type 1 - Large")
        gtkComboBoxAppendText(sw2, "Model 14: Realistic model Type 2 - Small")
        gtkComboBoxAppendText(sw2, "Model 15: Realistic model Type 2 - Medium")
        gtkComboBoxAppendText(sw2, "Model 16: Realistic model Type 2 - Large")
        gtkComboBoxAppendText(sw2, "Model 17: Realistic model Type 3 - Small")
        gtkComboBoxAppendText(sw2, "Model 18: Realistic model Type 3 - Medium")
        gtkComboBoxAppendText(sw2, "Model 19: Realistic model Type 3 - Large")
        sw3 <- gtkComboBoxNewText(show = TRUE)
        gtkComboBoxAppendText(sw3, "Interval of current electrodes, 1")
        gtkComboBoxAppendText(sw3, "Interval of current electrodes, 2")
        gtkComboBoxAppendText(sw3, "Interval of current electrodes, 3")
        gtkComboBoxAppendText(sw3, "Interval of current electrodes, 4")
        gtkComboBoxAppendText(sw3, "Interval of current electrodes, 5")
        gtkComboBoxAppendText(sw3, "Interval of current electrodes, 6")
        gtkComboBoxAppendText(sw3, "Interval of current electrodes, 7")
        gtkComboBoxAppendText(sw3, "Interval of current electrodes, 8")
        gtkComboBoxAppendText(sw3, "Interval of current electrodes, 9")
        gtkComboBoxAppendText(sw3, "Interval of current electrodes, 10")
        gtkComboBoxAppendText(sw3, "Interval of current electrodes, 11")
        gtkComboBoxAppendText(sw3, "Interval of current electrodes, 12")
        gtkComboBoxAppendText(sw3, "Interval of current electrodes, 13")
        gtkComboBoxAppendText(sw3, "Interval of current electrodes, 14")
        gtkComboBoxAppendText(sw3, "Interval of current electrodes, 15")
        sw4 <- gtkComboBoxNewText(show = TRUE)
        gtkComboBoxAppendText(sw4, "w/ estimation of subcutaneous fat & resolution enhacement")
        gtkComboBoxAppendText(sw4, "with estimation of subcutaneous fat area")
        gtkComboBoxAppendText(sw4, "w/o estimation of subcutaneous fat & resolution enhancement")
        gtkComboBoxAppendText(sw4, "without estimation of subcutaneous fat area")
        sw5 <- gtkComboBoxNewText(show = TRUE)
        gtkComboBoxAppendText(sw5, "Measurement with parallel pattern")
        gtkComboBoxAppendText(sw5, "Measurement with bypass pattern")
        vbox <- gtkVBoxNew()
        vbox$packStart(button1, expand = TRUE)
        vbox$packStart(button2, expand = TRUE)
        vbox$packStart(sw1, expand = FALSE)
        vbox$packStart(sw4, expand = FALSE)
        vbox$packStart(button3, expand = TRUE)
        vbox$packStart(button4, expand = TRUE)
        vbox$packStart(sw2, expand = FALSE)
        vbox$packStart(sw3, expand = FALSE)
        vbox$packStart(sw5, expand = FALSE)
        gSignalConnect(obj = button1, signal = "clicked", f = prep)
        gSignalConnect(obj = button2, signal = "clicked", f = recon)
        gSignalConnect(obj = button3, signal = "clicked", f = postProcess)
        gSignalConnect(obj = button4, signal = "clicked", f = sim)
        gtkComboBoxSetActive(sw1, 0)
        gtkComboBoxSetActive(sw2, 0)
        gtkComboBoxSetActive(sw3, 0)
        gtkComboBoxSetActive(sw4, 0)
        gtkComboBoxSetActive(sw5, 0)
        win$add(vbox)
        win$show()
    }
}
# END 

# Code takes less than 1 second to run
start.time <- proc.time()

# Reading the csv file with NYC Measles Data
read.ymdc <- function(file) {
  # Read the CSV file
  data <- read.csv(file, header = TRUE)
  data$Date <- as.Date(data$date.yyyy.mm.dd, format = "%Y-%m-%d")
  
  # Return the modified data
  return(data)
}


nyc <- read.ymdc("NYCWeeklyMeasles.csv") # Data
nyc <- nyc[, c("Date", "measles")]

colnames(nyc) <- c('Date','Cases') # Renaming columns


# Ensure the date column is in Date format
nyc$Cases[is.na(nyc$Cases)] <- mean(nyc$Cases, na.rm = TRUE) # get rid of missing data by using the mean
str(nyc) # confirm the type of data

nyc$Cases<-sqrt(nyc$Cases) # take the square root of the data to clearly see what is happening


# Function to plot the time series (with different colors for specific time ranges)
timeplot <- function(data, x, y, ma = NULL, add = FALSE, time_ranges=NULL, colors = NULL,...) {
  # Check if the specified columns exist in the data frame
  if (!(x %in% colnames(data)) | !(y %in% colnames(data))) {
    stop("The specified columns do not exist in the data frame.")
  }
  
  # Extract the data for plotting
  x <- data[[x]]
  y <- data[[y]]
  
  # Create the initial plot or add to it to existing plot if add is TRUE
  if (!add) {
    plot(x, y, type = "n",
         ylim = range(y, na.rm = TRUE), 
         xlim = range(x, na.rm = TRUE),
         ann = FALSE, bty = 'L', 
         xaxs = "i", las = 1, ...)
    lines(x, y, col = 'grey', ...)
  } else {
    lines(x, y, col = 'grey', ...)
  }
  
  # Add a moving average line if 'ma' is specified
  if (!is.null(ma)) {
    # Calculate the moving average
    y_ma <- filter(y, rep(1/ma, ma), sides = 2)
    lines(x, y_ma, col = 'grey', lwd = 2, ...)
  }
  
  if (!is.null(time_ranges) && !is.null(colors)) {
    for (i in seq_along(time_ranges)) {
      range <- time_ranges[[i]]
      color <- colors[i]
      
      # Subset data for the current time range
      in_range <- x >= as.Date(range[1]) & x <= as.Date(range[2])
      lines(x[in_range], y_ma[in_range], col = color, ...)
    }
  } else {
    # Plot the entire series in grey if no ranges are provided
    lines(x, y, col = 'grey', ...)
  }
}

# Function to plot periodogram
periodogram <- function(df, # data frame: date and (weekly) cases
                        xlim=c(0,10), # max 10 year period by default
                        trange, # whole time series by default
                        add=FALSE, color, # make new plot by default
                        ... ) {
  if (!is.na(trange[1])) { # consider only the specified time range
    df <- subset(df, Date >= min(trange) & Date <= max(trange))
  }
  with(df,{
    s <- spectrum(Cases, plot=FALSE)
    s$per <- 1/(52*s$freq) # period in years (assume data are weekly)
    if (!add) { # start a new plot without plotting data
      plot(s$per,s$spec,typ="n",bty="L",ann=FALSE,las=1,
           xaxs="i",yaxs="i", xlim=xlim, col=color,...)
    }
    lines( s$per, s$spec , col=color,lwd=2)
  })
}


#Multipanel plot for NYC Time Series Case

times<-list(c('1890-10-11','1945-12-29'),
            c('1946-1-5','1963-12-27'),
            c('1964-1-3','1982-12-30'))

colors<-c('hotpink','dodgerblue','seagreen2')

par(mfrow = c(2, 2), mar = c(5, 5, 4, 2))
timeplot(nyc, "Date", "Cases", ma = 6,time_ranges = times,
         colors = colors)
title(main = "Time Series of Cases", xlab = "Date", 
      ylab = "Number of Cases")
mtext('(a)',side=3,line=2,adj=0)

periodogram(nyc,trange=times[[1]],xlim=c(0,8),color=colors[[1]])
title(main = "Periodogram 1890-1945", xlab = "Years", 
      ylab = "Power Spectrum")
mtext('(b)',side=3,line=2,adj=0)

periodogram(nyc,trange=times[[2]],xlim=c(0,8),color=colors[[2]])
title(main = "Periodogram 1946-1963", xlab = "Years", 
      ylab = "Power Spectrum")
mtext('(c)',side=3,line=2,adj=0)

periodogram(nyc,trange=times[[3]],xlim=c(0,8),color=colors[[3]])
title(main = "Periodogram 1964-1980", xlab = "Years", 
      ylab = "Power Spectrum")
mtext('(d)',side=3,line=2,adj=0)

end.time <- proc.time() - start.time

print(end.time)

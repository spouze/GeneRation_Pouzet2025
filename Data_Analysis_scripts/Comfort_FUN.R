play = function(file=3){
  if (file==3){system("afplay /Users/sylvain/Music/snap.mp3")}
  if (file==1){system("afplay /Users/sylvain/Music/demisnap.mp3")}
  if (file==2){system("afplay /Users/sylvain/Music/demisnap2.mp3")}
  }

topologies = c("HIGHCO", "RANDOM", "SCALEF")
topocolors = c("red3", "green4", "blue4")

colrs = c("green4", "blue4", "red3")
# colrs[match("SCALEF", topologies)]

resetplot = function(){
  dev.off()
  dev.new()
  cex_mtext <<- 1
  par(mgp=c(3, 1, 0))
}

assign_category <- function(value, thresholds=effect_class_tresh, names=effect_class) {
  index <- findInterval(value, thresholds, rightmost.closed = T)
  return(names[index + 1])  # Adding 1 since R indexing starts at 1
}

check_folder_exists = function(folder_to_check){ # Check if the folder exists
  if (!dir.exists(folder_to_check)) { 
    dir.create(folder_to_check, recursive = TRUE)
    cat("Folder created at:", folder_to_check)
  } else{
    cat("Folder exists.")
  }
}

draw_box_check = function(){
  if (draw_box){
    box("figure", col="forestgreen")
    box(which="plot", col="red")
    box("outer", col="blue")
    }
}

restore_mgp = function(){
  par(mgp=c(3, 1, 0))
}

cmtoin=2.54 # 1 in = 2.54 cm

make_darker_color = function(color="#B3E4FF", amount=30){
  return(rgb(pmax(col2rgb(color) - amount, 0)[1,], pmax(col2rgb(color) - amount, 0)[2,], pmax(col2rgb(color) - amount, 0)[3,], maxColorValue = 255))
}

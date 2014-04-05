set terminal pngcairo 
set output 'serial_stream.png'
plot 'ss' using 1:2

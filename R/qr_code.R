# create a qr code
library(qrcode)
QR <- qr_code("https://benoitlepage.github.io/mediation_workshop/")
plot(QR)

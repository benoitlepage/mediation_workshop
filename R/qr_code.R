# create a qr code
library(qrcode)
QR <- qr_code("https://benoitlepage.github.io/mediation_workshop/")
plot(QR)


QR <- qr_code("https://drive.google.com/drive/folders/19EjXOO0rS0N4IdJ2XJo7Cb_j-lNpbNjs?usp=drive_link")
plot(QR)

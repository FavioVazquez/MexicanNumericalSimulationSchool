require(data.table)
require(rCharts)

df = fread("prueba1")type: 'logarithmic'

n1 <- hPlot(V2 ~ V1, data = df, type = "line")
#n1$xAxis(tickFormat = "#!function (x) { return Math.round(100*Math.pow(10,x))/100;}!#")
n1$yAxis(list(type = "logarithmic"))
n1

data(wuhan.hp)
whhp.hgwr <- hgwr(
  formula = Price ~ L(Fee + d.Water + d.Commercial +
    d.PrimarySchool + d.ShoppingMall + d.Kindergarten) +
    BuildingArea + (Floor.High | group),
  data = wuhan.hp,
  bw = "CV", kernel = "bisquared"
)

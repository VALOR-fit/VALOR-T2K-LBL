

SKexp new
Skobs new = realData

fparm new
PhysParamList [ntoys]

surface = th2d(nbins,nbins)

for (toy in ntoys):
  PhysParamList [itoy]->RandomizeAllSystematics()

for param1 in nbins:
  fparm->SetInitValue(param2)
  for param2 in nbins:
    fparm->Update(param2)
    likelihood = 0
   
    for itoy in ntoys:
      for ip in params:
        if not fparm->Fixed(i):
          PhysParamList[itoy].SetParam(ip)
          
      skexp->Generate(PhysParamList[itoy])
      chisq = skobs.Statistic(skexp)
      likelihood += exp(-chisq/2)
    #end toys
    
    final_likelihood = likelihood/ntoys
    chisq = -2*log(final_likelihood)
    
    surface[param1][param2].Set(param1,param2,chisq)

surface->Draw()

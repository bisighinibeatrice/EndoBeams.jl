Kᶜ¹¹ = [1,2,3]
hcat(vcat(Kᶜ¹¹, Kᶜ¹¹, Kᶜ¹¹, Kᶜ¹¹), vcat(Kᶜ¹¹, Kᶜ¹¹, Kᶜ¹¹, Kᶜ¹¹), vcat(Kᶜ¹¹, Kᶜ¹¹, Kᶜ¹¹,Kᶜ¹¹), vcat(Kᶜ¹¹, Kᶜ¹¹, Kᶜ¹¹, Kᶜ¹¹))

aux1 = 1:12
aux2 = 1:12
aux3 = 1:12
aux4 = 1:12

B = Mat224(
    aux1[1], aux3[1], 
    aux1[2], aux3[2], 
    aux1[3], aux3[3], 
    aux1[4], aux3[4], 
    aux1[5], aux3[5], 
    aux1[6], aux3[6],
    aux1[7], aux3[7], 
    aux1[8], aux3[8], 
    aux1[9], aux3[9], 
    aux1[10], aux3[10], 
    aux1[11], aux3[11], 
    aux1[12], aux3[12], 
    aux2[1], aux4[1], 
    aux2[2], aux4[2], 
    aux2[3], aux4[3], 
    aux2[4], aux4[4], 
    aux2[5], aux4[5], 
    aux2[6], aux4[6],
    aux2[7], aux4[7], 
    aux2[8], aux4[8], 
    aux2[9], aux4[9], 
    aux2[10], aux4[10], 
    aux2[11], aux4[11], 
    aux2[12], aux4[12])
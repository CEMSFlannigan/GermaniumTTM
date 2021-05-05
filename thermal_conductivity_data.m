temp = [3.25522483628520;3.47214698852705;3.43747584627629;3.63999854417698;3.82774140297810;4.03968181780061;4.26334242741817;4.47115049025887;4.76817668766587;5.10252003573688;5.45262978910836;5.85918537879986;6.37692432117432;6.90212166474628;7.47076260592536;8.11074067996867;8.75328745584954;9.44663786391612;10.2563618661359;11.1025208518409;12.0738016640040;12.9283353270309;14.1554240509142;15.1775807959248;16.1379601269318;17.4887294405275;18.9338952110449;20.3546379919692;20.3763075688763;21.9483673542492;23.6835252509822;25.0229961678234;26.5834767418582;28.7814425892344;31.1729845886038;32.6709623583556;33.6743499646640;34.8527984807871;37.5988262151247;40.5613006064517;43.9572768790624;46.7263481529062;47.2058262069529;50.0154023612886;53.1840117103856;56.5532747109717;59.9196349463222;62.3007349063629;65.7954152333920;68.4904303130104;72.5667608993591;77.0220829561716;81.4564263854966;84.4596228529386;89.8099913483561;95.1554814906713;100.196666792408;106.110397251621;111.553002925591;119.479590882361;127.509688441509;136.078484637662;144.660404456672;153.775621217632;161.267603107055;172.724372972635;184.329399344243;196.290901197206;209.932023942106;224.036962533947;239.089002286149;255.154810626589;271.318722242229;287.375122925286;304.337657310860;316.919469856505;340.656631806427;363.542937889388;388.648359889710;415.461113412006;420.052514629844;435.313647098817;468.042245044220;465.691254661238;481.669845434720;517.374726436794;517.984881674452;525.610211260345;558.598199916829;560.719503838829;570.601483466966;606.983260909786;609.579673375641;612.505860393730;647.160440382166;633.865457671153;664.886725641830;666.331164585943;714.729095733088;714.184492787646;742.231416094684;776.066925820504;794.580863955831;830.972835399316;883.840019106482;902.393741505416;950.927605100255;974.623094787426;1024.55850469163;1049.39871494491;1048.21739295227;1126.73937157256;1152.45398127063;1222.58619021507;1244.27134417718;1297.58434439506]; % K
therm_cond = [5.48742986871950;5.82139761613591;5.33405431815535;6.56594393722273;7.12256467036059;7.74388731368373;8.43404482016335;9.35104156802441;9.60825399795769;10.9020235774113;12.1851914566601;12.5017477939683;13.9300367776566;14.9614595758355;15.8680615361789;16.6575967360339;17.0172157274473;17.4765850996737;17.9198116160787;18.0296891635646;17.7871147919189;17.7039721362785;17.3522463311909;17.5469373806832;16.6866569278385;15.6003712955295;14.7621891860333;12.8607492679483;13.7714507114915;12.4693825466443;12.1978488802127;11.8634039739959;11.1251094494906;10.2995083122574;9.36184127488836;9.23827130579537;9.16104969667764;8.45986322380165;7.82291663399675;7.22601537251915;6.64199466984496;6.73906118821407;6.10509830786314;5.66265389109822;5.22822136568916;4.83078717390452;4.46407806219688;4.09361265220288;4.09817499366603;3.75940973898207;3.48824250126911;3.14739686904579;3.10150254177253;2.85502639022662;2.64147218449203;2.44725133613827;2.19372891273134;2.16104305753691;1.98987942273613;1.82710759147984;1.67349749119499;1.53840239663447;1.35758269240921;1.33545438282510;1.22992915521289;1.13708188894446;1.05166282444012;0.968566084265713;0.897408384566510;0.829993904190619;0.768577492528682;0.708253131564571;0.654911215142853;0.595193662161783;0.581401355112633;0.540600697264556;0.501934484929332;0.465358631306669;0.423138167827379;0.432586274738551;0.395673415249252;0.371150593309109;0.383789917224623;0.346275223943519;0.322434054806362;0.337124706052360;0.302476769057813;0.274031396695637;0.333844578795692;0.304221029147165;0.268377972778231;0.305108062956848;0.246955754347281;0.269671456085779;0.229206216269165;0.291949635286686;0.210973825569578;0.251899408065014;0.214061553308651;0.242154998365601;0.190786964062023;0.218390095323522;0.192180765736523;0.179033692829942;0.188421455029358;0.172164532452948;0.190827881902263;0.165777664034064;0.190287567643831;0.174857082816847;0.163476804708571;0.175448782643765;0.163524582333243;0.179731614814938;0.165732491930665;0.173465794794252]; % W/cm/K

log_temp = log(temp)/log(10);
log_therm_cond = log(therm_cond)/log(10);

log_therm_cond_fit = fit(log_temp,log_therm_cond,'smoothingspline','SmoothingParam',0.9996275803392751);

scatter(log_temp,log_therm_cond);
hold on;
plot(log_temp,log_therm_cond_fit(log_temp));

p1Ki =     0.01207;
p2Ki =     0.04158;
p3Ki =    -0.00384;
p4Ki =     -0.1541;
p5Ki =      0.1266;
p6Ki =      0.0513;
p7Ki =      -1.026;
p8Ki =      0.3893;
log_therm_cond_fit = @(LT) p1Ki*LT^7 + p2Ki*LT^6 + p3Ki*LT^5 + p4Ki*LT^4 + p5Ki*LT^3 + p6Ki*LT^2 + p7Ki*LT^1 + p8Ki; % (LT - 1.976)/0.7905
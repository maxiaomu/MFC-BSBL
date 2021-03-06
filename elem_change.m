function change_data = elem_change(h1)
[ch,hh] = size(h1);
jiange1 = 1:2:ch;
ch_ban = ch/2;
jiange2 = 1:1:ch_ban;
m1 = h1(jiange1);
m2 = h1(jiange1+1);
m3 = max(m1,m2);
output = zeros(ch,hh);
output(jiange2*2,:) = m3(jiange2,:);
output(jiange2*2-1,:) = m3(jiange2,:);
change_data = output;

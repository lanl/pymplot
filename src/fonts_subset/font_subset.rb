
fonts=[\
'Arial','ArialBold','ArialItalic', \
'Courier\ New','Courier\ NewBold','Courier\ NewItalic', \
'Courier\ Prime','Courier\ PrimeBold','Courier\ PrimeItalic', \
'Helvetica','HelveticaBold','HelveticaItalic', \
'Consolas','ConsolasBold','ConsolasItalic', \
'Times\ New\ Roman','Times\ New\ RomanBold','Times\ New\ RomanItalic', \
'IBMPlexSans','IBMPlexSansBold','IBMPlexSansItalic' \
'CMU\ Serif','CMU\ SerifBold','CMU\ SerifItalic', \
]

require "fileutils"
fonts.each do |i|
	system("FONTFORGE_LANGUAGE=ff ./extract_subset ../fonts/#{i}.ttf ./#{i}.ttf &")
end 

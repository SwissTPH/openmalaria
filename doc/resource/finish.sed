# Devilish cunning sed script for wrapping Markdown output as HTML
# (from http://www.alleged.org.uk/2005/marky/)

1   {
    h
    i\
<!DOCTYPE html PUBLIC "-//W3C//DTD XHTML 1.0 Strict//EN"\
    "http://www.w3.org/TR/xhtml1/DTD/xhtml1-strict.dtd">\
<html xml:lang="en-GB" lang="en-GB">\
<head>
    s/h1>/title>/g
    rresource/metadata.inc
    a\
</head>\
<body>
    rresource/header.inc
}

2   {
    x
    p
    x
}

$   {
    rresource/footer.inc
    a\
</body>\
</html>
}

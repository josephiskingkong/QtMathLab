QT       += core gui

greaterThan(QT_MAJOR_VERSION, 4): QT += widgets printsupport

CONFIG += c++17

SOURCES += \
    main.cpp \
    mainwindow.cpp \
    muParser.cpp \
    muParserBase.cpp \
    muParserBytecode.cpp \
    muParserCallback.cpp \
    muParserDLL.cpp \
    muParserError.cpp \
    muParserInt.cpp \
    muParserTest.cpp \
    muParserTokenReader.cpp \
    navbar.cpp \
    qcustomplot.cpp \
    secondwindow.cpp

HEADERS += \
    mainwindow.h \
    muParser.h \
    muParserBase.h \
    muParserBytecode.h \
    muParserCallback.h \
    muParserDLL.h \
    muParserDef.h \
    muParserError.h \
    muParserFixes.h \
    muParserInt.h \
    muParserTemplateMagic.h \
    muParserTest.h \
    muParserToken.h \
    muParserTokenReader.h \
    navbar.h \
    qcustomplot.h \
    secondwindow.h

FORMS += \
    mainwindow.ui \
    navbar.ui \
    secondwindow.ui

# Default rules for deployment.
qnx: target.path = /tmp/$${TARGET}/bin
else: unix:!android: target.path = /opt/$${TARGET}/bin
!isEmpty(target.path): INSTALLS += target

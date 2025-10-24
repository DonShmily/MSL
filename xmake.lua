set_project("MSL")
set_version("1.0.0")
set_xmakever("3.0.3")
set_warnings("all")
set_allowedplats("windows", "linux", "macosx", "mingw")

add_rules("mode.debug", "mode.release")
set_config("plat", "mingw")
if is_plat("mingw") then
    set_config("sdk", "C:/Programing/msys64/ucrt64")
    set_toolchains("gcc")

elseif is_plat("linux") then
    set_toolchains("gcc")
end

set_languages("c++20")

add_requires("eigen3", {system = true})

add_includedirs("include")

if is_plat("mingw") then
        set_targetdir("$(projectdir)/bin/mingw",{ bindir = "bin", libdir = "lib" })
    elseif is_plat("windows") then
        set_targetdir("$(projectdir)/bin/windows",{ bindir = "bin", libdir = "lib" })
    elseif is_plat("linux") then
        set_targetdir("$(projectdir)/bin/linux",{ bindir = "bin", libdir = "lib" })
    end

add_cxflags("-fPIC")

includes("tests")
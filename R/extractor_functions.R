## Extract hexbin cell IDs on the fly for any supported object
setGeneric(".extract_cID", function(object) standardGeneric(".extract_cID"))
setMethod(".extract_cID", "ANY", function(object) {
    # compute and return cell assignments
    get_hexbin_info(object)$cID
})

## Extract hexbin matrix on the fly for any supported object
setGeneric(".extract_hexbin", function(object) standardGeneric(".extract_hexbin"))
setMethod(".extract_hexbin", "ANY", function(object) {
    # compute and return hexbin matrix
    get_hexbin_info(object)$hexbin.matrix
})

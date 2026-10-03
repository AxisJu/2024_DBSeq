def compute(inputs):
    """Compute using supplied in-memory input mappings, helper functions and options."""
    import sys
    from PIL import Image
    left, right = [inputs['left_right_1'][str(Image)].convert('RGB') for p in sys.argv[1:3]]
    out = Image.new('RGB', (left.width + right.width, max(left.height, right.height)), 'white')
    out.paste(left, (0, 0))
    out.paste(right, (left.width, 0))
    out.save(sys.argv[3])
    return locals()
